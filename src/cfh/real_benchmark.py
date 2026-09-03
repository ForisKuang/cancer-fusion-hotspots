"""Prototype live cBioPortal-to-domain-retention benchmark pipeline."""

from __future__ import annotations

import csv
import json
import math
import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.algorithms.registry import list_algorithms
from cfh.genes.registry import GeneConfig, load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.domain_source import ProteinDomain
from cfh.mapping.feature_mapper import map_event
from cfh.mapping.genome_nexus_source import GenomeNexusClient, resolve_domains
from cfh.mapping.transcript_source import resolve_breakpoint_protein_position
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.normalization.event_normalizer import normalize
from cfh.orchestrator.run import run_algorithms
from cfh.stats.breakpoint_tests import build_frame_domain_contingency_table


class RealBenchmarkError(RuntimeError):
    """Base class for expected, user-actionable benchmark failures."""


class RealBenchmarkInputError(RealBenchmarkError):
    """Raised when the requested gene/study cannot define a valid run."""


class RealBenchmarkNetworkError(RealBenchmarkError):
    """Raised when a required public data service cannot be reached."""


@dataclass
class RealBenchmarkRun:
    gene_symbol: str
    study_id: str
    molecular_profile_id: str
    retrieved_at: datetime
    raw_structural_variant_count: int
    events: list[FusionEvent]
    features: list[FusionFeature]
    rows: list[dict]
    results: list[AlgorithmResult]
    summary: dict
    warnings: list[str]
    endpoints: list[str] = field(default_factory=list)
    reference: dict | None = None


class _ResolvedDomainSource:
    """Expose already-resolved Genome Nexus domains to ``map_event``."""

    def __init__(self, domains: list[ProteinDomain]):
        self.domains = domains

    def fetch(self, _accession: str) -> list[ProteinDomain]:
        return self.domains


def _is_target_protein_fusion(event: FusionEvent, target_gene: str) -> bool:
    genes = {str(event.Site1_gene or "").upper(), str(event.Site2_gene or "").upper()}
    return (
        target_gene.upper() in genes
        and event.Is_protein_fusion is True
        and "protein fusion" in str(event.Event_info or "").lower()
    )


def _target_breakpoint(row: dict, target_gene: str) -> int:
    target = target_gene.upper()
    if str(row.get("Site1_Hugo_Symbol") or "").upper() == target:
        value = row.get("Site1_Position")
    elif str(row.get("Site2_Hugo_Symbol") or "").upper() == target:
        value = row.get("Site2_Position")
    else:  # pragma: no cover - guarded by _is_target_protein_fusion
        raise ValueError(f"row does not contain target gene {target_gene}")
    if value is None or pd.isna(value):
        raise ValueError(f"{target_gene} fusion has no genomic breakpoint")
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{target_gene} fusion has invalid genomic breakpoint {value!r}") from exc


def _target_role(event: FusionEvent, target_gene: str) -> str:
    target = target_gene.upper()
    if str(event.Five_prime_gene or "").upper() == target:
        return "five_prime"
    if str(event.Three_prime_gene or "").upper() == target:
        return "three_prime"
    raise ValueError(
        f"could not determine 5'/3' role for {target_gene} in {event.Event_id}; "
        f"Event_Info={event.Event_info!r}"
    )


def _partner(event: FusionEvent, target_gene: str) -> str:
    target = target_gene.upper()
    if str(event.Site1_gene or "").upper() == target:
        return event.Site2_gene or "unknown"
    if str(event.Site2_gene or "").upper() == target:
        return event.Site1_gene or "unknown"
    return "unknown"


def _source_frame_status(text: object) -> str | None:
    """Read an explicit frame call without applying production normalization rules."""
    value = str(text or "").strip().lower()
    if not value:
        return None
    if value in {"na", "n/a", "unknown"}:
        return None
    if "out-of-frame" in value or "out of frame" in value:
        return "out-of-frame"
    if "in-frame" in value or "in frame" in value:
        return "in-frame"
    return None


def _load_benchmark_config(gene_symbol: str) -> GeneConfig:
    try:
        config = load_gene_config(gene_symbol)
    except FileNotFoundError as exc:
        raise RealBenchmarkInputError(
            f"Unknown gene {gene_symbol!r}. Run `cfh list-genes` to see configured genes."
        ) from exc
    if config.gene_symbol is None or config.entrez_gene_id is None:
        raise RealBenchmarkInputError(
            f"Gene {gene_symbol!r} needs gene_symbol and entrez_gene_id in its config."
        )
    if not config.key_domains:
        raise RealBenchmarkInputError(
            f"Gene {gene_symbol!r} has no key domain configured for domain-retention analysis."
        )
    return config


def _unavailable_domain_result(
    message: str,
    events: list[FusionEvent],
    features: list[FusionFeature],
    config: GeneConfig,
    n_permutations: int,
) -> AlgorithmResult:
    return AlgorithmResult(
        Algorithm="domain_retention",
        Algorithm_version="0.1.0",
        Parameters={"seed": 42, "n_permutations": n_permutations},
        Summary={
            "fisher_odds_ratio": None,
            "fisher_p_value": None,
            "permutation_empirical_p_value": None,
            "observed_in_frame_retention_rate": None,
        },
        Tables={
            "frame_domain_contingency_table": build_frame_domain_contingency_table(
                events, features, config
            ),
            "permutation_null_retention_rates": [],
        },
        Warnings=[message],
    )


def analyze_structural_variant_calls(
    calls: list[dict],
    gene_symbol: str,
    study_id: str,
    *,
    molecular_profile_id: str | None = None,
    genome_nexus_client: GenomeNexusClient | None = None,
    n_permutations: int = 1_000,
    algorithm_names: list[str] | None = None,
) -> RealBenchmarkRun:
    """Normalize and analyze already-fetched cBioPortal SV API objects."""
    if n_permutations <= 0:
        raise RealBenchmarkInputError("n_permutations must be positive")
    config = _load_benchmark_config(gene_symbol)
    profile_id = molecular_profile_id or f"{study_id}_structural_variants"
    client = genome_nexus_client or GenomeNexusClient()

    raw = cbioportal_api.structural_variants_to_dataframe(calls)
    normalized = normalize(raw, None, study_id)
    selected = [
        (row.to_dict(), event)
        for (_, row), event in zip(raw.iterrows(), normalized, strict=True)
        if _is_target_protein_fusion(event, config.gene_symbol)
    ]

    warnings: list[str] = []
    if selected:
        try:
            domains = resolve_domains(
                config.gene_symbol,
                config.protein_id,
                genome_nexus_client=client,
            )
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus/domain lookup failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
    else:
        domains = []
        warnings.append(
            f"No Protein Fusion records for {config.gene_symbol} were returned by "
            f"{profile_id}. Verify the gene and study ID, or try a study with SV data."
        )
    domain_source = _ResolvedDomainSource(domains)
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    rows: list[dict] = []
    target_key = config.key_domains[0].key or config.key_domains[0].name

    for row, event in selected:
        try:
            role = _target_role(event, config.gene_symbol)
            breakpoint = _target_breakpoint(row, config.gene_symbol)
            mapping = resolve_breakpoint_protein_position(
                None,
                config,
                breakpoint_genomic=breakpoint,
                genome_nexus_client=client,
            )
            if mapping.breakpoint_protein_position is None:
                raise ValueError("Genome Nexus returned no protein position")
            feature = map_event(
                event,
                config,
                role=role,
                junction_position_aa=mapping.breakpoint_protein_position,
                domain_source=domain_source,
            ).model_copy(update={"Breakpoint_exon": mapping.breakpoint_exon})
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus breakpoint mapping failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
        except Exception as exc:
            warnings.append(
                f"Skipped {event.Event_id} ({event.Fusion_name or 'unnamed fusion'}): "
                f"{type(exc).__name__}: {exc}"
            )
            continue
        events.append(event)
        features.append(feature)
        rows.append(
            {
                "event_id": event.Event_id,
                "sample_id": event.Sample_id,
                "patient_id": event.Patient_id,
                "fusion_name": event.Fusion_name,
                "partner_gene": _partner(event, config.gene_symbol),
                "frame_status": event.Frame_status,
                "target_role": role,
                "breakpoint_genomic": breakpoint,
                "breakpoint_exon": mapping.breakpoint_exon,
                "breakpoint_protein_position": mapping.breakpoint_protein_position,
                "is_intronic_breakpoint": mapping.is_intronic_breakpoint,
                "domain_status": (feature.Domain_retention_flags or {}).get(target_key, "unknown"),
                "retained_domains": "; ".join(feature.Retained_domains or []),
                "lost_domains": "; ".join(feature.Lost_domains or []),
                "disrupted_domains": "; ".join(feature.Disrupted_domains or []),
                "source_annotation_text": " | ".join(
                    value
                    for value in (
                        str(row.get("Annotation") or ""),
                        str(row.get("Event_Info") or ""),
                    )
                    if value
                ),
                "source_site2_effect_on_frame": row.get("Site2_Effect_On_Frame"),
            }
        )

    in_frame_mapped = any(
        event.Frame_status == "in-frame"
        and (feature.Domain_retention_flags or {}).get(target_key)
        in {"retained", "lost", "disrupted"}
        for event, feature in zip(events, features, strict=True)
    )
    if not in_frame_mapped:
        message = (
            "Domain-retention statistics are unavailable because no mapped in-frame "
            "protein-fusion record has a known domain state. Verify the gene/study IDs "
            "and source frame annotations."
        )
        warnings.append(message)
        domain_result = _unavailable_domain_result(
            message, events, features, config, n_permutations
        )
    else:
        try:
            domain_result = run_algorithms(
                ["domain_retention"],
                events,
                features,
                config,
                {"domain_retention": {
                    "seed": 42,
                    "n_permutations": n_permutations,
                    "genome_nexus_client": client,
                }},
            )[0]
        except requests.RequestException as exc:
            raise RealBenchmarkNetworkError(
                f"Genome Nexus permutation mapping failed for {config.gene_symbol}: {exc}. "
                "Check network access and https://www.genomenexus.org availability, then retry."
            ) from exc
        except ValueError as exc:
            message = f"Domain-retention statistics are unavailable: {exc}"
            warnings.append(message)
            domain_result = _unavailable_domain_result(
                message, events, features, config, n_permutations
            )
    selected_events = [event for _, event in selected]
    requested_algorithms = algorithm_names or ["domain_retention", "frequency"]
    other_algorithms = [name for name in requested_algorithms if name != "domain_retention"]
    other_results = run_algorithms(
        other_algorithms,
        selected_events,
        features,
        config,
        {
            "frequency": {"dedup_by_patient": False},
            "cutpoint_detection": {
                "n_permutations": n_permutations,
                "genome_nexus_client": client,
            },
        },
    )
    results_by_name = {result.Algorithm: result for result in [domain_result, *other_results]}
    results = [results_by_name[name] for name in requested_algorithms if name in results_by_name]
    frequency_result = results_by_name.get("frequency")
    if frequency_result is None:
        frequency_result = FrequencyAnalysis().run(
            selected_events, features, config, {"dedup_by_patient": False}
        )
    partner_counts = frequency_result.Tables["Partner_gene_counts"]
    in_frame_count = sum(event.Frame_status == "in-frame" for event in selected_events)
    retained_count = sum(row["domain_status"] == "retained" for row in rows)
    in_frame_retained_count = sum(
        row["frame_status"] == "in-frame" and row["domain_status"] == "retained"
        for row in rows
    )
    domain_definition = next(
        (
            domain
            for domain in domains
            if domain.accession == config.key_domains[0].accession
            or domain.name == config.key_domains[0].accession
        ),
        None,
    )
    total = len(selected_events)
    summary = {
        "raw_structural_variant_count": len(calls),
        "total_fusions": total,
        "mapped_fusions": len(features),
        "skipped_fusions": total - len(features),
        "in_frame_count": in_frame_count,
        "in_frame_percent": 100 * in_frame_count / total if total else 0.0,
        "kinase_retained_count": retained_count,
        "kinase_retained_percent": 100 * retained_count / total if total else 0.0,
        "in_frame_kinase_retained_count": in_frame_retained_count,
        "fisher_odds_ratio": domain_result.Summary["fisher_odds_ratio"],
        "fisher_p_value": domain_result.Summary["fisher_p_value"],
        "permutation_p_value": domain_result.Summary["permutation_empirical_p_value"],
        "frame_domain_contingency_table": domain_result.Tables[
            "frame_domain_contingency_table"
        ],
        "partner_counts": partner_counts,
        "domain_accession": config.key_domains[0].accession,
        "domain_start_aa": domain_definition.start_aa if domain_definition else None,
        "domain_end_aa": domain_definition.end_aa if domain_definition else None,
    }
    return RealBenchmarkRun(
        gene_symbol=config.gene_symbol,
        study_id=study_id,
        molecular_profile_id=profile_id,
        retrieved_at=datetime.now(timezone.utc),
        raw_structural_variant_count=len(calls),
        events=events,
        features=features,
        rows=rows,
        results=results,
        summary=summary,
        warnings=warnings,
        endpoints=[
            "https://www.cbioportal.org/api/structural-variant/fetch",
            "https://www.genomenexus.org/ensembl/canonical-transcript/hgnc/"
            f"{config.gene_symbol}",
        ],
        reference=(
            config.benchmark_reference.model_dump() if config.benchmark_reference else None
        ),
    )


def run_real_benchmark(
    gene_symbol: str,
    study_id: str,
    *,
    n_permutations: int = 1_000,
    algorithm_names: list[str] | None = None,
) -> RealBenchmarkRun:
    """Fetch and analyze a gene's structural variants from cBioPortal."""
    config = _load_benchmark_config(gene_symbol)
    profile_id = f"{study_id}_structural_variants"
    try:
        calls = cbioportal_api.fetch_structural_variants(
            [config.entrez_gene_id],
            [profile_id],
        )
    except requests.RequestException as exc:
        raise RealBenchmarkNetworkError(
            f"cBioPortal request failed for gene {config.gene_symbol} and profile "
            f"{profile_id}: {exc}. Check the study ID, network access, and "
            "https://www.cbioportal.org availability, then retry."
        ) from exc
    return analyze_structural_variant_calls(
        calls,
        gene_symbol,
        study_id,
        molecular_profile_id=profile_id,
        n_permutations=n_permutations,
        algorithm_names=algorithm_names,
    )


def run_analysis(
    gene_symbol: str,
    study_id: str,
    *,
    n_permutations: int = 1_000,
) -> RealBenchmarkRun:
    """Run every registered plugin against the shared live mapped input snapshot."""
    return run_real_benchmark(
        gene_symbol,
        study_id,
        n_permutations=n_permutations,
        algorithm_names=list_algorithms(),
    )


def _json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _format_stat(value: float | int | None) -> str:
    return "unavailable" if value is None else f"{value:.6g}"


def markdown_summary(run: RealBenchmarkRun) -> str:
    """Render a concise, checked-in-friendly benchmark report."""
    summary = run.summary
    domain = summary["domain_accession"] or "configured domain"
    partners = ", ".join(
        f"{row['Partner_gene']} ({row['Event_count']})" for row in summary["partner_counts"]
    )
    table = summary["frame_domain_contingency_table"]
    lines = [
        f"# {run.gene_symbol} real-data fusion benchmark: {run.study_id}",
        "",
        "Retrieved from public cBioPortal and Genome Nexus on "
        f"{run.retrieved_at.date().isoformat()}.",
        "",
        "## Results",
        "",
        f"- Structural variants returned for {run.gene_symbol}: "
        f"{summary['raw_structural_variant_count']}",
        f"- Protein-fusion records found: {summary['total_fusions']}",
        f"- Protein-fusion records mapped: {summary['mapped_fusions']}",
        f"- Malformed/unmappable fusion records skipped: {summary['skipped_fusions']}",
        f"- In-frame: {summary['in_frame_count']}/{summary['total_fusions']} "
        f"({summary['in_frame_percent']:.1f}%)",
        f"- {domain} ({summary['domain_start_aa']}-{summary['domain_end_aa']} aa) retained: "
        f"{summary['kinase_retained_count']}/{summary['total_fusions']} "
        f"({summary['kinase_retained_percent']:.1f}%)",
        f"- In-frame and {domain}-retained: "
        f"{summary['in_frame_kinase_retained_count']}/{summary['in_frame_count']}",
        f"- Fisher exact test (one-sided): odds ratio "
        f"{_format_stat(summary['fisher_odds_ratio'])}, "
        f"p={_format_stat(summary['fisher_p_value'])}",
        "- Breakpoint-permutation empirical p-value: "
        f"{_format_stat(summary['permutation_p_value'])}",
        f"- Contingency table `[[retained/in-frame, retained/other], "
        f"[not-retained/in-frame, not-retained/other]]`: `{table}`",
        "",
        "## Method",
        "",
        f"The cBioPortal `{run.molecular_profile_id}` structural-variant profile was "
        "queried by the configured Entrez gene ID. Records explicitly annotated as "
        "`Protein Fusion` were adapted to the production SV schema and normalized; "
        "the raw `site2EffectOnFrame=NA` values were therefore resolved from "
        "`Event_Info`, not copied into `FusionEvent.Frame_status`.",
        "",
        f"{run.gene_symbol} genomic breakpoints were mapped against the Genome Nexus "
        f"canonical transcript, and retention was classified against its returned {domain} "
        "coordinates. Counts are event-level with no patient deduplication. The "
        "Fisher comparison's `other` column combines out-of-frame and unknown-frame "
        "events, as pre-specified by the domain-retention algorithm.",
        "",
        "## Partners",
        "",
        partners or "None",
        "",
    ]
    if run.warnings:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in run.warnings)
        lines.append("")
    if run.reference:
        citation = run.reference["citation"]
        lines.extend(
            [
                "## Reference comparison",
                "",
                f"| Metric | {citation} | This run |",
                "|---|---:|---:|",
                f"| In-frame | {run.reference['in_frame_percent']:.1f}% | "
                f"{summary['in_frame_percent']:.1f}% |",
                f"| Domain retained | {run.reference['domain_retained_percent']:.1f}% | "
                f"{summary['kinase_retained_percent']:.1f}% |",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            "",
        ]
    )
    if not summary["total_fusions"]:
        lines.append("No fusion records were available, so comparison is not possible.")
    elif run.gene_symbol.upper() == "BRAF" and run.study_id == "msk_impact_50k_2026":
        lines.extend(
            [
                "This does **not** reproduce the Zehir et al. (PMC5461196) report of "
                "33/33 BRAF fusions being in-frame with the kinase domain retained: "
                f"this live successor cohort has {summary['in_frame_count']}/"
                f"{summary['total_fusions']} in-frame and "
                f"{summary['in_frame_kinase_retained_count']}/"
                f"{summary['in_frame_count']} in-frame fusions retaining {domain}.",
                "",
                "`msk_impact_50k_2026` is a newer successor cohort, not the paper's "
                "original `msk_impact_2017` cohort. This is therefore replication in a "
                "related cohort, not a reanalysis of the paper's original 33 cases.",
            ]
        )
    else:
        lines.append("These values describe the live study named above.")
    lines.append("")
    return "\n".join(lines)


def _git_sha() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _run_id(run: RealBenchmarkRun) -> str:
    gene = re.sub(r"[^a-z0-9]+", "-", run.gene_symbol.lower()).strip("-")
    study = re.sub(r"[^a-z0-9]+", "-", run.study_id.lower()).strip("-")
    timestamp = run.retrieved_at.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{gene}_{study}_{timestamp}"


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    names = fieldnames or (list(rows[0]) if rows else [])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n")
        if names:
            writer.writeheader()
            writer.writerows(rows)


def _discrepancies(run: RealBenchmarkRun) -> list[dict]:
    discrepancies = []
    for row in run.rows:
        common = {
            "event_id": row.get("event_id"),
            "partner_gene": row.get("partner_gene"),
            "frame_status": row.get("frame_status"),
            "retained_domains": row.get("retained_domains", ""),
            "lost_domains": row.get("lost_domains", ""),
            "disrupted_domains": row.get("disrupted_domains", ""),
            "source_annotation_text": row.get("source_annotation_text", ""),
            "source_site2_effect_on_frame": row.get("source_site2_effect_on_frame"),
        }
        if run.reference and (
            row.get("frame_status") != "in-frame" or row.get("domain_status") != "retained"
        ):
            discrepancies.append({"discrepancy_type": "reference_discrepancy", **common})
        source_calls = {
            status
            for status in (
                _source_frame_status(row.get("source_site2_effect_on_frame")),
                _source_frame_status(row.get("source_annotation_text")),
            )
            if status in {"in-frame", "out-of-frame"}
        }
        if any(source_call != row.get("frame_status") for source_call in source_calls):
            discrepancies.append({"discrepancy_type": "source_vs_derived_qa_mismatch", **common})
    return discrepancies


def _domain_track_svg(run: RealBenchmarkRun, outlier_ids: set[str]) -> str:
    """Render mapped breakpoints and configured-domain boundaries without gene literals."""
    start = run.summary.get("domain_start_aa") or 0
    end = run.summary.get("domain_end_aa") or 0
    positions = [
        row["breakpoint_protein_position"]
        for row in run.rows
        if row["breakpoint_protein_position"]
    ]
    maximum = max([end, *positions, 1])
    scale = 800 / maximum
    dots = []
    for index, row in enumerate(run.rows):
        position = row["breakpoint_protein_position"]
        if position is None:
            continue
        color = "#d62728" if row["event_id"] in outlier_ids else "#2878b5"
        y = 100 + (index % 5) * 7
        dots.append(f'<circle cx="{60 + position * scale:.1f}" cy="{y}" r="3" fill="{color}"/>')
    return "\n".join([
        '<svg xmlns="http://www.w3.org/2000/svg" width="920" height="180" viewBox="0 0 920 180">',
        '<rect width="920" height="180" fill="white"/>',
        f'<text x="60" y="28" font-family="sans-serif" font-size="16">'
        f'{run.gene_symbol} domain-retention track</text>',
        '<line x1="60" y1="90" x2="860" y2="90" stroke="#444" stroke-width="4"/>',
        f'<rect x="{60 + start * scale:.1f}" y="76" '
        f'width="{max(2, (end-start)*scale):.1f}" height="28" '
        'fill="#62b36f" opacity="0.55"/>',
        *dots,
        '<circle cx="60" cy="157" r="4" fill="#2878b5"/>'
        '<text x="70" y="162" font-family="sans-serif" font-size="12">'
        'matches reference pattern</text>',
        '<circle cx="270" cy="157" r="4" fill="#d62728"/>'
        '<text x="280" y="162" font-family="sans-serif" font-size="12">'
        'reference discrepancy</text>',
        '</svg>',
    ])


def _comparison_svg(run: RealBenchmarkRun) -> str:
    reference = run.reference or {}
    metrics = [
        ("In-frame", reference.get("in_frame_percent", 0), run.summary["in_frame_percent"]),
        (
            "Domain retained",
            reference.get("domain_retained_percent", 0),
            run.summary["kinase_retained_percent"],
        ),
    ]
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="620" height="210" '
        'viewBox="0 0 620 210">',
        '<rect width="620" height="210" fill="white"/>',
        f'<text x="20" y="25" font-family="sans-serif" font-size="16">'
        f'Reference vs {run.study_id}</text>',
    ]
    for index, (label, ref_value, run_value) in enumerate(metrics):
        y = 55 + index * 70
        elements.extend([
            f'<text x="20" y="{y}" font-family="sans-serif" font-size="12">{label}</text>',
            f'<rect x="140" y="{y-14}" width="{ref_value*3.8:.1f}" height="16" fill="#999"/>',
            f'<text x="530" y="{y}" font-family="sans-serif" font-size="12">'
            f'reference {ref_value:.1f}%</text>',
            f'<rect x="140" y="{y+10}" width="{run_value*3.8:.1f}" height="16" fill="#2878b5"/>',
            f'<text x="530" y="{y+24}" font-family="sans-serif" font-size="12">'
            f'run {run_value:.1f}%</text>',
        ]
        )
    elements.append("</svg>")
    return "\n".join(elements)


def write_outputs(
    run: RealBenchmarkRun,
    output_dir: str | Path,
    *,
    output_stem: str | None = None,
    cli_args: list[str] | None = None,
    run_id: str | None = None,
) -> dict[str, Path]:
    """Write a complete, provenance-bearing run directory."""
    del output_stem  # retained as a compatibility-only keyword for older callers
    destination = Path(output_dir) / (run_id or _run_id(run))
    visualization_dir = destination / "visualizations"
    visualization_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = destination / "results.tsv"
    json_path = destination / "results.json"
    markdown_path = destination / "report.md"

    _write_tsv(tsv_path, run.rows)

    algorithm_results = [result.model_dump(mode="json") for result in run.results]
    for result in algorithm_results:
        null_rates = (result.get("Tables") or {}).get("permutation_null_retention_rates")
        if null_rates is not None:
            result["Tables"]["permutation_null_retention_rates"] = {
                "omitted_from_artifact": True,
                "count": len(null_rates),
            }
    payload = _json_safe(
        {
            "gene_symbol": run.gene_symbol,
            "study_id": run.study_id,
            "molecular_profile_id": run.molecular_profile_id,
            "retrieved_at": run.retrieved_at.isoformat(),
            "summary": run.summary,
            "warnings": run.warnings,
            "events": run.rows,
            "algorithm_results": algorithm_results,
        }
    )
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    markdown_path.write_text(markdown_summary(run))
    discrepancies = _discrepancies(run)
    outliers_path = destination / "outliers.tsv"
    _write_tsv(outliers_path, discrepancies, [
        "discrepancy_type", "event_id", "partner_gene", "frame_status",
        "retained_domains", "lost_domains", "disrupted_domains",
        "source_annotation_text", "source_site2_effect_on_frame",
    ])
    reference_ids = {
        row["event_id"] for row in discrepancies
        if row["discrepancy_type"] == "reference_discrepancy"
    }
    domain_svg = visualization_dir / "domain_retention_outliers.svg"
    comparison_svg = visualization_dir / "reference_comparison.svg"
    domain_svg.write_text(_domain_track_svg(run, reference_ids) + "\n")
    comparison_svg.write_text(_comparison_svg(run) + "\n")
    manifest_path = destination / "manifest.json"
    manifest_path.write_text(json.dumps({
        "gene": run.gene_symbol,
        "study_id": run.study_id,
        "endpoints_used": run.endpoints,
        "git_sha": _git_sha(),
        "cli_args": cli_args or [],
        "timestamp": run.retrieved_at.isoformat(),
    }, indent=2) + "\n")
    return {
        "run_directory": destination,
        "manifest": manifest_path,
        "tsv": tsv_path,
        "json": json_path,
        "markdown": markdown_path,
        "outliers": outliers_path,
        "domain_svg": domain_svg,
        "comparison_svg": comparison_svg,
    }

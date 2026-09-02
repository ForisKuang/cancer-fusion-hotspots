"""Prototype live cBioPortal-to-domain-retention benchmark pipeline."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.domain_source import ProteinDomain
from cfh.mapping.feature_mapper import map_event
from cfh.mapping.genome_nexus_source import GenomeNexusClient, resolve_domains
from cfh.mapping.transcript_source import resolve_breakpoint_protein_position
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.normalization.event_normalizer import normalize


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
    if value is None:
        raise ValueError(f"{target_gene} fusion has no genomic breakpoint")
    return int(value)


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


def analyze_structural_variant_calls(
    calls: list[dict],
    gene_symbol: str,
    study_id: str,
    *,
    molecular_profile_id: str | None = None,
    genome_nexus_client: GenomeNexusClient | None = None,
    n_permutations: int = 1_000,
) -> RealBenchmarkRun:
    """Normalize and analyze already-fetched cBioPortal SV API objects."""
    config = load_gene_config(gene_symbol)
    if config.gene_symbol is None or config.entrez_gene_id is None:
        raise ValueError(f"{gene_symbol} config lacks gene_symbol or entrez_gene_id")
    profile_id = molecular_profile_id or f"{study_id}_structural_variants"
    client = genome_nexus_client or GenomeNexusClient()

    raw = cbioportal_api.structural_variants_to_dataframe(calls)
    normalized = normalize(raw, None, study_id)
    selected = [
        (row.to_dict(), event)
        for (_, row), event in zip(raw.iterrows(), normalized, strict=True)
        if _is_target_protein_fusion(event, config.gene_symbol)
    ]

    domains = resolve_domains(
        config.gene_symbol,
        config.protein_id,
        genome_nexus_client=client,
    )
    domain_source = _ResolvedDomainSource(domains)
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    rows: list[dict] = []
    target_key = config.key_domains[0].key or config.key_domains[0].name

    for row, event in selected:
        role = _target_role(event, config.gene_symbol)
        breakpoint = _target_breakpoint(row, config.gene_symbol)
        mapping = resolve_breakpoint_protein_position(
            None,
            config,
            breakpoint_genomic=breakpoint,
            genome_nexus_client=client,
        )
        if mapping.breakpoint_protein_position is None:
            raise ValueError(f"Genome Nexus did not map a protein position for {event.Event_id}")
        feature = map_event(
            event,
            config,
            role=role,
            junction_position_aa=mapping.breakpoint_protein_position,
            domain_source=domain_source,
        ).model_copy(update={"Breakpoint_exon": mapping.breakpoint_exon})
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
            }
        )

    domain_result = DomainRetentionAlgorithm().run(
        events,
        features,
        config,
        {
            "seed": 42,
            "n_permutations": n_permutations,
            "genome_nexus_client": client,
        },
    )
    frequency_result = FrequencyAnalysis().run(
        events, features, config, {"dedup_by_patient": False}
    )
    results = [domain_result, frequency_result]
    partner_counts = frequency_result.Tables["Partner_gene_counts"]
    in_frame_count = sum(event.Frame_status == "in-frame" for event in events)
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
    total = len(events)
    summary = {
        "raw_structural_variant_count": len(calls),
        "total_fusions": total,
        "mapped_fusions": len(features),
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
    )


def run_real_benchmark(
    gene_symbol: str,
    study_id: str,
    *,
    n_permutations: int = 1_000,
) -> RealBenchmarkRun:
    """Fetch and analyze a gene's structural variants from cBioPortal."""
    config = load_gene_config(gene_symbol)
    if config.entrez_gene_id is None:
        raise ValueError(f"{gene_symbol} config lacks entrez_gene_id")
    profile_id = f"{study_id}_structural_variants"
    calls = cbioportal_api.fetch_structural_variants(
        [config.entrez_gene_id],
        [profile_id],
    )
    return analyze_structural_variant_calls(
        calls,
        gene_symbol,
        study_id,
        molecular_profile_id=profile_id,
        n_permutations=n_permutations,
    )


def _json_safe(value: object) -> object:
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


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
        f"- Protein-fusion records analyzed: {summary['total_fusions']}",
        f"- In-frame: {summary['in_frame_count']}/{summary['total_fusions']} "
        f"({summary['in_frame_percent']:.1f}%)",
        f"- {domain} ({summary['domain_start_aa']}-{summary['domain_end_aa']} aa) retained: "
        f"{summary['kinase_retained_count']}/{summary['total_fusions']} "
        f"({summary['kinase_retained_percent']:.1f}%)",
        f"- In-frame and {domain}-retained: "
        f"{summary['in_frame_kinase_retained_count']}/{summary['in_frame_count']}",
        f"- Fisher exact test (one-sided): odds ratio "
        f"{summary['fisher_odds_ratio']:.6g}, p={summary['fisher_p_value']:.6g}",
        f"- Breakpoint-permutation empirical p-value: {summary['permutation_p_value']:.6g}",
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
        "## Interpretation",
        "",
    ]
    if run.gene_symbol.upper() == "BRAF" and run.study_id == "msk_impact_50k_2026":
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


def write_outputs(
    run: RealBenchmarkRun,
    output_dir: str | Path,
    *,
    output_stem: str | None = None,
) -> dict[str, Path]:
    """Write event-level TSV, structured JSON, and a Markdown summary."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    stem = output_stem or f"{run.gene_symbol.lower()}_{run.study_id}_benchmark"
    tsv_path = destination / f"{stem}.tsv"
    json_path = destination / f"{stem}.json"
    markdown_path = destination / f"{stem}.md"

    fieldnames = list(run.rows[0]) if run.rows else []
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        if fieldnames:
            writer.writeheader()
            writer.writerows(run.rows)

    algorithm_results = [result.model_dump(mode="json") for result in run.results]
    for result in algorithm_results:
        null_rates = (result.get("Tables") or {}).get("permutation_null_retention_rates")
        if null_rates is not None:
            result["Tables"]["permutation_null_retention_rates"] = {
                "omitted_from_artifact": True,
                "count": len(null_rates),
            }
    payload = _json_safe({
        "gene_symbol": run.gene_symbol,
        "study_id": run.study_id,
        "molecular_profile_id": run.molecular_profile_id,
        "retrieved_at": run.retrieved_at.isoformat(),
        "summary": run.summary,
        "events": run.rows,
        "algorithm_results": algorithm_results,
    })
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    markdown_path.write_text(markdown_summary(run))
    return {"tsv": tsv_path, "json": json_path, "markdown": markdown_path}

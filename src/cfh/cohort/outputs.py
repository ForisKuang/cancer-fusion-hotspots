"""Write ``cfh cohort-scan`` output artifacts to ``runs/<run_id>/cohort_scan/``.

Two kinds of artifact are produced:

* One consolidated summary (TSV + JSON + Markdown + PDF) covering every
  scanned gene -- recurrence, in-frame%, domain-retention%, FDR-adjusted
  p-values, and composite score -- sorted by significance. The summary also
  calls out an "honorable mentions" tier (see :func:`build_honorable_mentions`):
  the top-ranked genes by raw Fisher p-value among genes that did NOT survive
  genome-wide FDR correction, so a human reviewer can see the real
  second-tier signal that a strict significant/not-significant split hides.
* Full detailed per-gene reports (the same ``results.tsv``/``results.json``/
  ``report.md``/``report.pdf``/``manifest.json`` that ``cfh analyze`` already
  produces, via :func:`cfh.real_benchmark.write_outputs`) for every
  hand-curated gene, every FDR-significant (q<0.05) gene, and every
  honorable-mention gene -- never for every scanned gene, so a run with
  hundreds of non-significant genes does not generate hundreds of PDFs.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from cfh.cohort.scan import (
    CohortScanResult,
    genes_needing_full_report,
    per_gene_min_q_value,
)
from cfh.real_benchmark import write_outputs
from cfh.reporting.manhattan import render_manhattan_svg
from cfh.reporting.pdf import render_cohort_summary_pdf

_MANHATTAN_SVG_FILENAME = "manhattan.svg"

DEFAULT_HONORABLE_MENTION_COUNT = 15
"""Default size of the "honorable mentions" / highly ranked non-FDR-significant
tier: the top N genes by raw Fisher p-value among genes that did NOT survive
genome-wide FDR correction. Configurable via ``write_cohort_scan_outputs``'s
``honorable_mention_count`` (and the ``cfh cohort-scan
--honorable-mention-count`` CLI flag)."""

_HONORABLE_MENTION_HEADING = (
    "Honorable mentions: highly ranked non-FDR-significant genes worth human review"
)

_HONORABLE_MENTION_NOTE = (
    "Did not survive genome-wide multiple-testing correction (FDR-adjusted "
    "q-value at or above the significance threshold), but ranks highly by "
    "raw p-value among the non-FDR-significant genes and may warrant "
    "targeted follow-up. This is NOT a claim of statistical significance."
)

_SUMMARY_FIELDNAMES = [
    "gene_symbol",
    "config_source",
    "status",
    "distinct_patient_count",
    "total_sv_count",
    "n_events_analyzed",
    "in_frame_percent",
    "domain_retention_percent",
    "fisher_p_value",
    "permutation_p_value",
    "min_fdr_adjusted_q_value",
    "fdr_significant",
    "top_composite_score",
    "top_composite_partner_gene",
    "error",
]


def _top_composite(run) -> tuple[float | None, str | None]:
    if run is None:
        return None, None
    for result in run.results:
        if result.Algorithm != "composite_score":
            continue
        ranking = (result.Tables or {}).get("composite_evidence_ranking") or []
        scored = [row for row in ranking if row.get("Composite_score") is not None]
        if not scored:
            return None, None
        top = max(scored, key=lambda row: row["Composite_score"])
        return top["Composite_score"], top.get("Partner_gene")
    return None, None


def build_summary_rows(result: CohortScanResult) -> list[dict]:
    """Build one row per scanned gene, sorted by significance (most
    significant FDR-adjusted q-value first; genes with no q-value sort
    last, tie-broken by descending recurrence).
    """
    min_q_by_gene = per_gene_min_q_value(result)
    significant = set(result.significant_genes)
    rows: list[dict] = []
    for outcome in result.gene_outcomes:
        summary = outcome.run.summary if outcome.run is not None else {}
        top_score, top_partner = _top_composite(outcome.run)
        q_value = min_q_by_gene.get(outcome.gene_symbol)
        rows.append(
            {
                "gene_symbol": outcome.gene_symbol,
                "config_source": outcome.config_source,
                "status": outcome.status,
                "distinct_patient_count": outcome.distinct_patient_count,
                "total_sv_count": outcome.total_sv_count,
                "n_events_analyzed": summary.get("total_fusions"),
                "in_frame_percent": summary.get("in_frame_percent"),
                "domain_retention_percent": summary.get("kinase_retained_percent"),
                "fisher_p_value": summary.get("fisher_p_value"),
                "permutation_p_value": summary.get("permutation_p_value"),
                "min_fdr_adjusted_q_value": q_value,
                "fdr_significant": outcome.gene_symbol in significant,
                "top_composite_score": top_score,
                "top_composite_partner_gene": top_partner,
                "error": outcome.error,
            }
        )

    def _sort_key(row: dict) -> tuple:
        q_value = row["min_fdr_adjusted_q_value"]
        return (
            q_value is None,
            q_value if q_value is not None else 1.0,
            -(row["distinct_patient_count"] or 0),
            row["gene_symbol"],
        )

    rows.sort(key=_sort_key)
    return rows


def build_honorable_mentions(
    rows: list[dict], *, limit: int = DEFAULT_HONORABLE_MENTION_COUNT
) -> list[dict]:
    """The "honorable mentions" / highly ranked non-FDR-significant tier: the
    top ``limit`` genes by raw Fisher p-value among genes that are NOT
    FDR-significant.

    This is a distinct, additive ranking -- never a claim that any of these
    genes IS significant, or even that it is close to significant (a run's
    544 genes can rank from a raw p-value of 0.0004 down to 0.9, and this
    tier surfaces the top of that non-significant ranking regardless of how
    far any individual gene actually sits from the threshold). A run can
    genuinely have only one FDR-significant gene (a strict, correct
    genome-wide multiple-testing result) while still having a real,
    biologically sensible second tier of genes that rank highly on raw
    p-value but did not survive correction; without this, that tier is
    invisible to a human reviewer skimming a binary significant/not-significant
    summary.
    """
    candidates = [
        row
        for row in rows
        if not row["fdr_significant"] and row.get("fisher_p_value") is not None
    ]
    candidates.sort(key=lambda row: (row["fisher_p_value"], row["gene_symbol"]))
    mentions = []
    for rank, row in enumerate(candidates[:limit], start=1):
        mentions.append(
            {
                "rank": rank,
                "gene_symbol": row["gene_symbol"],
                "fisher_p_value": row["fisher_p_value"],
                "min_fdr_adjusted_q_value": row["min_fdr_adjusted_q_value"],
                "n_events_analyzed": row["n_events_analyzed"],
                "in_frame_percent": row["in_frame_percent"],
                "domain_retention_percent": row["domain_retention_percent"],
                "note": _HONORABLE_MENTION_NOTE,
            }
        )
    return mentions


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _format_cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


_HONORABLE_MENTION_FIELDNAMES = [
    "rank",
    "gene_symbol",
    "fisher_p_value",
    "min_fdr_adjusted_q_value",
    "n_events_analyzed",
    "in_frame_percent",
    "domain_retention_percent",
]


def _write_summary_markdown(
    result: CohortScanResult,
    rows: list[dict],
    honorable_mentions: list[dict],
    path: Path,
) -> None:
    lines = [
        f"# Genome-wide fusion-hotspot cohort scan: {result.study_id}",
        "",
        f"- Total genes with any structural-variant record in the cohort: "
        f"{result.total_genes_before_gating}",
        f"- Genes passing the >= {result.min_distinct_patients}-distinct-patient "
        f"recurrence gate: {result.genes_after_gating}",
        f"- Curated gene configs used: {result.curated_gene_count}",
        f"- Auto-generated gene configs used: {result.auto_config_gene_count}",
        f"- Genes gated in but unresolvable (no Genome Nexus canonical transcript): "
        f"{result.unresolved_gene_count}",
        f"- FDR-significant genes (q < 0.05) after Benjamini-Hochberg correction across "
        f"all {len(rows)} scanned genes: {len(result.significant_genes)}",
        "",
    ]
    if result.warnings:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in result.warnings)
        lines.append("")
    if honorable_mentions:
        lines.extend(
            [
                f"## {_HONORABLE_MENTION_HEADING}",
                "",
                f"The following {len(honorable_mentions)} gene(s) **did not survive genome-wide "
                "multiple-testing correction** (FDR-adjusted q-value at or above the "
                f"q={result.significance_level:g} significance threshold), but rank highest by "
                "raw Fisher p-value among the non-FDR-significant genes and may warrant targeted "
                "follow-up. This section is **not** a claim of statistical significance -- see "
                "the FDR-significant genes above for that.",
                "",
            ]
        )
        lines.append("| " + " | ".join(_HONORABLE_MENTION_FIELDNAMES) + " |")
        lines.append("|" + "---|" * len(_HONORABLE_MENTION_FIELDNAMES))
        for mention in honorable_mentions:
            cells = " | ".join(
                _format_cell(mention[field]) for field in _HONORABLE_MENTION_FIELDNAMES
            )
            lines.append(f"| {cells} |")
        lines.append("")
    lines.extend(
        [
            "## Genome-wide summary plot",
            "",
            f"One point per scanned gene with an FDR-adjusted q-value: x-axis is genes "
            "ranked by composite evidence score (descending), y-axis is "
            "-log10(FDR-adjusted q-value), with a dashed line at the "
            f"q={result.significance_level:g} significance threshold.",
            "",
            f"![Genome-wide fusion-hotspot summary plot]({_MANHATTAN_SVG_FILENAME})",
            "",
        ]
    )
    lines.extend(["## Scanned genes (sorted by significance)", ""])
    lines.append("| " + " | ".join(_SUMMARY_FIELDNAMES) + " |")
    lines.append("|" + "---|" * len(_SUMMARY_FIELDNAMES))
    for row in rows:
        cells = " | ".join(_format_cell(row[field]) for field in _SUMMARY_FIELDNAMES)
        lines.append(f"| {cells} |")
    lines.append("")
    path.write_text("\n".join(lines))


def _write_summary_pdf(
    result: CohortScanResult, rows: list[dict], honorable_mentions: list[dict], path: Path
) -> None:
    header = [
        "Gene",
        "Source",
        "Status",
        "Patients",
        "In-frame%",
        "Domain-ret.%",
        "Best q",
        "Sig.",
        "Composite",
    ]
    table_rows = [header]
    for row in rows:
        table_rows.append(
            [
                row["gene_symbol"],
                row["config_source"],
                row["status"],
                _format_cell(row["distinct_patient_count"]),
                _format_cell(row["in_frame_percent"]),
                _format_cell(row["domain_retention_percent"]),
                _format_cell(row["min_fdr_adjusted_q_value"]),
                "yes" if row["fdr_significant"] else "no",
                _format_cell(row["top_composite_score"]),
            ]
        )

    extra_tables = []
    if honorable_mentions:
        mention_header = [
            "Rank",
            "Gene",
            "Fisher p",
            "Best q",
            "N events",
            "In-frame%",
            "Domain-ret.%",
        ]
        mention_rows = [mention_header] + [
            [
                _format_cell(mention["rank"]),
                mention["gene_symbol"],
                _format_cell(mention["fisher_p_value"]),
                _format_cell(mention["min_fdr_adjusted_q_value"]),
                _format_cell(mention["n_events_analyzed"]),
                _format_cell(mention["in_frame_percent"]),
                _format_cell(mention["domain_retention_percent"]),
            ]
            for mention in honorable_mentions
        ]
        extra_tables.append(
            {
                "heading": _HONORABLE_MENTION_HEADING,
                "note": (
                    "Did NOT survive genome-wide FDR correction -- ranked by raw Fisher "
                    "p-value among non-significant genes. Not a claim of significance."
                ),
                "rows": mention_rows,
            }
        )

    render_cohort_summary_pdf(
        path,
        title=f"Genome-wide fusion-hotspot cohort scan: {result.study_id}",
        subtitle=(
            f"{result.genes_after_gating} of {result.total_genes_before_gating} genes passed the "
            f">= {result.min_distinct_patients}-patient recurrence gate; "
            f"{len(result.significant_genes)} FDR-significant (q<0.05)."
        ),
        notes=[f"Generated {datetime.now(timezone.utc).isoformat()}"],
        rows=table_rows,
        extra_tables=extra_tables,
        figures_dir=path.parent,
    )


def write_cohort_scan_outputs(
    result: CohortScanResult,
    output_dir: str | Path,
    *,
    run_id: str | None = None,
    pdf: bool = True,
    honorable_mention_count: int = DEFAULT_HONORABLE_MENTION_COUNT,
) -> dict[str, Path]:
    """Write the consolidated summary and full per-gene reports for a
    completed cohort scan to ``<output_dir>/<run_id>/cohort_scan/``."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = run_id or f"cohort-scan_{result.study_id}_{timestamp}"
    destination = Path(output_dir) / run_id / "cohort_scan"
    destination.mkdir(parents=True, exist_ok=True)

    rows = build_summary_rows(result)
    honorable_mentions = build_honorable_mentions(rows, limit=honorable_mention_count)

    curated_genes = {
        outcome.gene_symbol
        for outcome in result.gene_outcomes
        if outcome.config_source == "curated"
    }
    honorable_mention_genes = {mention["gene_symbol"] for mention in honorable_mentions}
    priority_genes = curated_genes | honorable_mention_genes

    manhattan_svg_path = destination / _MANHATTAN_SVG_FILENAME
    manhattan_svg_path.write_text(
        render_manhattan_svg(
            rows, significance_level=result.significance_level, priority_genes=priority_genes
        )
        + "\n"
    )

    tsv_path = destination / "summary.tsv"
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=_SUMMARY_FIELDNAMES, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    json_path = destination / "summary.json"
    payload = _json_safe(
        {
            "study_id": result.study_id,
            "min_distinct_patients": result.min_distinct_patients,
            "total_genes_before_gating": result.total_genes_before_gating,
            "genes_after_gating": result.genes_after_gating,
            "curated_gene_count": result.curated_gene_count,
            "auto_config_gene_count": result.auto_config_gene_count,
            "unresolved_gene_count": result.unresolved_gene_count,
            "significant_genes": result.significant_genes,
            "warnings": result.warnings,
            "genes": rows,
            "honorable_mentions": honorable_mentions,
        }
    )
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    markdown_path = destination / "summary.md"
    _write_summary_markdown(result, rows, honorable_mentions, markdown_path)

    paths = {
        "run_directory": destination,
        "summary_tsv": tsv_path,
        "summary_json": json_path,
        "summary_markdown": markdown_path,
        "manhattan_svg": manhattan_svg_path,
    }

    if pdf:
        pdf_path = destination / "summary.pdf"
        _write_summary_pdf(result, rows, honorable_mentions, pdf_path)
        paths["summary_pdf"] = pdf_path

    gene_reports_dir = destination / "gene_reports"
    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    full_report_genes = genes_needing_full_report(result, honorable_mention_genes)
    gene_report_paths: dict[str, dict[str, Path]] = {}
    for gene_symbol in full_report_genes:
        outcome = outcomes_by_gene.get(gene_symbol)
        if outcome is None or outcome.run is None:
            continue
        gene_report_paths[gene_symbol] = write_outputs(
            outcome.run,
            gene_reports_dir,
            run_id=gene_symbol.lower(),
            pdf=pdf,
            cli_args=["cohort-scan", result.study_id],
        )
    paths["gene_reports"] = gene_report_paths

    return paths

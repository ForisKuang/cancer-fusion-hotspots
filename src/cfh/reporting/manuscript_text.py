"""Deterministic, ``summary.json``-derived text templates for the cross-gene
manuscript-style synthesis report (``paper.md``/``paper.pdf``).

Every sentence produced here is computed directly from numeric/list fields
already present in a cohort scan's ``summary.json`` payload -- the same dict
shape written by ``cfh.cohort.outputs.write_cohort_scan_outputs`` (the
``genes`` rows from ``build_summary_rows`` and the ``honorable_mentions``
list from ``build_honorable_mentions``), plus a small number of additive
top-level fields (``significance_level``, ``generated_at``). Nothing here
calls an LLM and nothing here is a fixed generic sentence independent of the
run's numbers: a needed value that is absent causes the corresponding
sentence/clause to be omitted or replaced with an explicit
"unavailable"/"unknown number of" statement, never an invented value --
exactly the discipline already established by
:func:`cfh.reporting.text.render_abstract`, which this module extends
rather than reimplements (the shared numeric formatting helpers are reused
directly from there).

Every function here is a pure function of its arguments, so the produced
text is byte-for-byte reproducible for a given ``summary.json``.
"""

from __future__ import annotations

from cfh.algorithms.registry import list_algorithms
from cfh.reporting.text import format_percent, format_stat, significance_clause


def _plural(count: int | None, singular: str = "", plural: str = "s") -> str:
    return singular if count == 1 else plural


def _count_display(value: int | None) -> str:
    return "an unknown number of" if value is None else str(value)


def _level_display(significance_level: float | None) -> str:
    return "the configured" if significance_level is None else f"{significance_level:g}"


def render_manuscript_title(payload: dict) -> str:
    """Render the manuscript's title from the scan's study id."""
    study = payload.get("study_id") or "the configured study"
    return f"Genome-wide fusion-hotspot analysis of {study}"


def render_manuscript_abstract(payload: dict) -> str:
    """Render the 3-6 sentence ABSTRACT section for one cohort scan's
    ``summary.json``-shaped payload."""
    study = payload.get("study_id") or "the configured study"
    total_before = payload.get("total_genes_before_gating")
    total_after = payload.get("genes_after_gating")
    min_patients = payload.get("min_distinct_patients")
    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    unresolved = payload.get("unresolved_gene_count")

    sentences: list[str] = []

    if total_before is not None and total_after is not None and min_patients is not None:
        sentences.append(
            f"This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of {study}: "
            f"{total_before} gene{_plural(total_before)} carried at least one structural-variant "
            f"record in the cohort, of which {total_after} passed the >= {min_patients}-distinct-"
            f"patient recurrence gate and were analyzed with the full registered algorithm suite "
            f"({_count_display(curated)} using hand-curated gene configs, {_count_display(auto)} "
            f"auto-configured, {_count_display(unresolved)} gated in but unresolvable)."
        )
    else:
        sentences.append(
            f"This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of {study}."
        )

    significance_level = payload.get("significance_level")
    level_display = _level_display(significance_level)
    significant_genes = payload.get("significant_genes") or []
    rows_by_gene = {row["gene_symbol"]: row for row in payload.get("genes") or []}

    if significant_genes:

        def _q_sort_key(gene: str) -> float:
            q_value = rows_by_gene.get(gene, {}).get("min_fdr_adjusted_q_value")
            return q_value if q_value is not None else 1.0

        clauses = []
        for gene in sorted(significant_genes, key=_q_sort_key):
            row = rows_by_gene.get(gene, {})
            n_events = row.get("n_events_analyzed")
            clauses.append(
                f"{gene} ({_count_display(n_events)} event{_plural(n_events)}, "
                f"{format_percent(row.get('in_frame_percent'))} in-frame, "
                f"{format_percent(row.get('domain_retention_percent'))} domain-retained, "
                f"Fisher p={format_stat(row.get('fisher_p_value'))}, "
                f"q={format_stat(row.get('min_fdr_adjusted_q_value'))})"
            )
        sentences.append(
            f"{len(significant_genes)} gene{_plural(len(significant_genes))} reached genome-wide "
            f"Benjamini-Hochberg FDR significance (q < {level_display}): "
            + "; ".join(clauses)
            + "."
        )
    else:
        sentences.append(
            "No gene reached genome-wide Benjamini-Hochberg FDR significance "
            f"(q < {level_display}) in this scan."
        )

    honorable_mentions = payload.get("honorable_mentions") or []
    if honorable_mentions:
        verb = "forms" if len(honorable_mentions) == 1 else "form"
        sentences.append(
            f"{len(honorable_mentions)} additional gene{_plural(len(honorable_mentions))} {verb} "
            "a highly ranked non-FDR-significant tier flagged for targeted follow-up (see "
            "Honorable mentions, below)."
        )
    else:
        sentences.append("No honorable-mentions tier was produced for this scan.")

    return " ".join(sentences)


def render_manuscript_methods(payload: dict) -> str:
    """Render the templated METHODS paragraph: data source, gating
    threshold, the registered algorithm suite (read programmatically from
    :func:`cfh.algorithms.registry.list_algorithms`, never hardcoded), and
    the statistical tests applied."""
    study = payload.get("study_id") or "the configured study"
    min_patients = payload.get("min_distinct_patients")
    total_before = payload.get("total_genes_before_gating")
    total_after = payload.get("genes_after_gating")
    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    unresolved = payload.get("unresolved_gene_count")

    sentences: list[str] = []

    if min_patients is not None and total_before is not None and total_after is not None:
        sentences.append(
            f"Structural-variant records were retrieved from the {study} cBioPortal study and "
            f"gated to genes with at least {min_patients} distinct patient(s) carrying a "
            f"structural-variant record ({total_after} of {total_before} genes passed the gate: "
            f"{_count_display(curated)} using hand-curated gene configs and {_count_display(auto)} "
            "auto-configured from Genome Nexus canonical-transcript/Pfam-domain data, with "
            f"{_count_display(unresolved)} gated-in gene(s) left unresolvable)."
        )
    else:
        sentences.append(
            f"Structural-variant records were retrieved from the {study} cBioPortal study."
        )

    algorithm_names = list_algorithms()
    if algorithm_names:
        sentences.append(
            "Each gated gene was analyzed with the full registered algorithm suite "
            f"({', '.join(algorithm_names)})."
        )
    else:
        sentences.append("No algorithms were registered for this run.")

    n_genes_scanned = len(payload.get("genes") or [])
    level_display = _level_display(payload.get("significance_level"))
    sentences.append(
        "Domain-retention and domain-disruption significance were assessed per gene with "
        "Fisher's exact test and a breakpoint-position permutation test; the resulting p-values "
        f"across all {n_genes_scanned} scanned gene{_plural(n_genes_scanned)} were jointly "
        f"corrected with Benjamini-Hochberg false-discovery-rate correction at q < {level_display}."
    )
    return " ".join(sentences)


def render_manhattan_caption(payload: dict) -> str:
    """Render the templated caption for the embedded genome-wide Manhattan
    figure: how many genes are plotted, the significance threshold, and how
    many genes sit above it -- computed with the same filter
    ``cfh.reporting.manhattan._plottable_points`` applies, so the caption
    never drifts from what the figure actually shows."""
    rows = payload.get("genes") or []
    plottable = [
        row
        for row in rows
        if row.get("min_fdr_adjusted_q_value") is not None
        and row.get("top_composite_score") is not None
    ]
    above = [row for row in plottable if row.get("fdr_significant")]
    level_display = _level_display(payload.get("significance_level"))
    return (
        f"Genome-wide summary of {len(plottable)} scanned gene{_plural(len(plottable))} with an "
        "FDR-adjusted q-value, ranked left-to-right by composite evidence score; the dashed line "
        f"marks the q={level_display} significance threshold, with {len(above)} "
        f"gene{_plural(len(above))} above it."
    )


def render_gene_highlight(row: dict, *, honorable_mention_note: str | None = None) -> str:
    """Render the 2-5 sentence templated highlight paragraph for one
    highlighted gene's summary row: gene name, event count, in-frame%,
    domain-retention%, and Fisher/FDR-adjusted p-values -- each in its own
    sentence, explicitly labeled with the statistic (raw p vs. BH-adjusted
    q) its significance verdict describes, since a gene's raw p-value can be
    significant while its genome-wide q-value is not -- following the same
    omission-over-invention discipline as
    :func:`cfh.reporting.text.render_abstract`. ``honorable_mention_note``,
    if given, is the gene's own already-established honorable-mentions note
    text (see :func:`cfh.cohort.outputs.build_honorable_mentions`), appended
    verbatim rather than re-derived.
    """
    gene = row.get("gene_symbol") or "This gene"
    n_events = row.get("n_events_analyzed")
    sentences: list[str] = []

    if n_events is not None:
        sentences.append(
            f"{gene} was analyzed across {n_events} fusion event{_plural(n_events)}, "
            f"{format_percent(row.get('in_frame_percent'))} in-frame and "
            f"{format_percent(row.get('domain_retention_percent'))} domain-retained."
        )
    else:
        sentences.append(f"Event counts are unavailable for {gene} in this scan.")

    fisher_p = row.get("fisher_p_value")
    fisher_display = format_stat(fisher_p)
    if fisher_display != "unavailable":
        sig = significance_clause(fisher_p)
        sentence = f"Domain-retention Fisher's exact test p={fisher_display}"
        sentence += f" (raw {sig})." if sig else "."
        sentences.append(sentence)

        q_display = format_stat(row.get("min_fdr_adjusted_q_value"))
        if q_display != "unavailable":
            # ``fdr_significant`` is the same precomputed BH-correction verdict
            # used for this gene's badge/tier elsewhere (see
            # ``cfh.cohort.outputs._gene_badges``) -- reusing it here (rather
            # than re-deriving a threshold check on q) guarantees this clause
            # can never drift from, or be confused with, the raw p-value's own
            # significance clause above. A gene's raw Fisher p-value can be
            # significant while its genome-wide BH-adjusted q-value is not
            # (that is the entire point of multiple-testing correction), so
            # these two verdicts are stated in separate sentences, each
            # explicitly labeled with the statistic it describes.
            fdr_significant = row.get("fdr_significant")
            if fdr_significant is True:
                q_clause = "reaches genome-wide FDR significance"
            elif fdr_significant is False:
                q_clause = "does not reach genome-wide FDR significance"
            else:
                q_clause = None
            q_sentence = f"Genome-wide BH-adjusted q={q_display}"
            q_sentence += f" ({q_clause})." if q_clause else "."
            sentences.append(q_sentence)
    else:
        sentences.append(
            f"No domain-retention statistical test could be computed for {gene} in this scan."
        )

    if honorable_mention_note:
        sentences.append(honorable_mention_note)

    return " ".join(sentences)


def render_discussion_bullets(payload: dict) -> list[str]:
    """Render the templated INTERPRETATION/DISCUSSION caveats. Each caveat
    that cites a number is gated on that number actually being present in
    ``payload`` (never invented); the final caveat is a standing
    interpretive statement -- true regardless of any specific number in this
    run -- not a claim derived from a numeric field.
    """
    bullets: list[str] = []

    n_genes_scanned = len(payload.get("genes") or [])
    level_display = _level_display(payload.get("significance_level"))
    bullets.append(
        "Cross-gene Benjamini-Hochberg FDR correction was applied jointly across all "
        f"{n_genes_scanned} scanned gene{_plural(n_genes_scanned)}' p-values. This reduces "
        "false-positive findings relative to testing each gene in isolation, but is a "
        f"conservative correction: a real per-gene effect can fail to reach the "
        f"q < {level_display} threshold once corrected across the full scanned gene set."
    )

    curated = payload.get("curated_gene_count")
    auto = payload.get("auto_config_gene_count")
    if auto:
        bullets.append(
            f"{_count_display(curated)} gene(s) used hand-curated gene configurations, while "
            f"{auto} gene(s) were auto-configured from Genome Nexus canonical-transcript/Pfam-"
            "domain data using a kinase/catalytic-keyword heuristic to select the tracked domain; "
            "auto-configured domains have not been manually verified the way hand-curated ones "
            "have."
        )

    generated_at = payload.get("generated_at")
    study = payload.get("study_id") or "the configured study"
    if generated_at:
        bullets.append(
            f"Data were retrieved live from the public {study} cBioPortal study as of "
            f"{generated_at}. Some cBioPortal cohorts (e.g. actively accruing clinical-"
            "sequencing panels such as MSK-IMPACT) are updated periodically, so exact counts "
            "could shift if this scan is re-run later against such a cohort."
        )

    bullets.append(
        "All findings in this manuscript are computational, hypothesis-generating candidate "
        "evidence from bioinformatic analysis of public cohort data, not validated clinical "
        "calls; a gene's presence here is not a therapeutic or diagnostic recommendation."
    )
    return bullets

"""Unit tests for the deterministic summary.json -> manuscript text templates
(:mod:`cfh.reporting.manuscript_text`).

Every assertion checks literal, byte-for-byte sentence text produced from a
synthetic ``summary.json``-shaped fixture, per the same never-LLM-generated,
never-a-fixed-generic-sentence discipline already enforced for
``cfh.reporting.text.render_abstract`` in ``tests/test_report_text.py``.
"""

from __future__ import annotations

from cfh.algorithms.registry import list_algorithms
from cfh.reporting.manuscript_text import (
    render_discussion_bullets,
    render_gene_highlight,
    render_manhattan_caption,
    render_manuscript_abstract,
    render_manuscript_methods,
    render_manuscript_title,
)

PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS = {
    "study_id": "msk_impact_50k_2026",
    "min_distinct_patients": 5,
    "total_genes_before_gating": 3919,
    "genes_after_gating": 544,
    "curated_gene_count": 4,
    "auto_config_gene_count": 520,
    "unresolved_gene_count": 20,
    "significant_genes": ["ETV6"],
    "significance_level": 0.05,
    "generated_at": "2026-09-04T14:36:53.295005+00:00",
    "genes": [
        {
            "gene_symbol": "ETV6",
            "config_source": "auto",
            "status": "ok",
            "n_events_analyzed": 90,
            "in_frame_percent": 71.11111111111111,
            "domain_retention_percent": 75.55555555555556,
            "fisher_p_value": 5.123255667283157e-06,
            "min_fdr_adjusted_q_value": 0.004334274294521551,
            "fdr_significant": True,
            "top_composite_score": 0.41577810494717626,
        },
        {
            "gene_symbol": "RET",
            "config_source": "curated",
            "status": "ok",
            "n_events_analyzed": 194,
            "in_frame_percent": 75.25773195876289,
            "domain_retention_percent": 92.26804123711341,
            "fisher_p_value": 0.00041966557652448966,
            "min_fdr_adjusted_q_value": 0.11976161828924137,
            "fdr_significant": False,
            "top_composite_score": 0.3,
        },
    ],
    "honorable_mentions": [
        {
            "rank": 1,
            "gene_symbol": "RET",
            "fisher_p_value": 0.00041966557652448966,
            "min_fdr_adjusted_q_value": 0.11976161828924137,
            "n_events_analyzed": 194,
            "in_frame_percent": 75.25773195876289,
            "domain_retention_percent": 92.26804123711341,
            "note": (
                "Did not survive genome-wide multiple-testing correction (FDR-adjusted "
                "q-value at or above the significance threshold), but ranks highly by raw "
                "p-value among the non-FDR-significant genes and may warrant targeted "
                "follow-up. This is NOT a claim of statistical significance."
            ),
        }
    ],
}

PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS = {
    "study_id": "some_other_study",
    "min_distinct_patients": 5,
    "total_genes_before_gating": 10,
    "genes_after_gating": 2,
    "curated_gene_count": 0,
    "auto_config_gene_count": 2,
    "unresolved_gene_count": 0,
    "significant_genes": [],
    "significance_level": 0.05,
    "generated_at": None,
    "genes": [
        {
            "gene_symbol": "FAKE1",
            "config_source": "auto",
            "status": "ok",
            "n_events_analyzed": 3,
            "in_frame_percent": 66.6,
            "domain_retention_percent": None,
            "fisher_p_value": None,
            "min_fdr_adjusted_q_value": None,
            "fdr_significant": False,
            "top_composite_score": None,
        }
    ],
    "honorable_mentions": [],
}


def test_render_manuscript_title_uses_study_id():
    assert (
        render_manuscript_title(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
        == "Genome-wide fusion-hotspot analysis of msk_impact_50k_2026"
    )


def test_manuscript_abstract_with_significant_and_honorable_mentions_matches_exact_text():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    assert abstract == (
        "This manuscript synthesizes a genome-wide fusion-hotspot cohort scan of "
        "msk_impact_50k_2026: 3919 genes carried at least one structural-variant record in "
        "the cohort, of which 544 passed the >= 5-distinct-patient recurrence gate and were "
        "analyzed with the full registered algorithm suite (4 using hand-curated gene "
        "configs, 520 auto-configured, 20 gated in but unresolvable). 1 gene reached "
        "genome-wide Benjamini-Hochberg FDR significance (q < 0.05): ETV6 (90 events, 71.1% "
        "in-frame, 75.6% domain-retained, Fisher p=5.12326e-06, q=0.00433427). 1 additional "
        "gene forms a highly ranked non-FDR-significant tier flagged for targeted follow-up "
        "(see Honorable mentions, below)."
    )


def test_manuscript_abstract_states_no_significant_genes_and_no_mentions_explicitly():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)

    assert (
        "No gene reached genome-wide Benjamini-Hochberg FDR significance (q < 0.05) in this "
        "scan." in abstract
    )
    assert "No honorable-mentions tier was produced for this scan." in abstract
    # Never invents a claim of significance or a mention count that isn't there.
    assert "0 additional" not in abstract


def test_manuscript_abstract_never_contains_overclaiming_honorable_mention_language():
    abstract = render_manuscript_abstract(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
    assert "near-significant" not in abstract.lower()
    assert "near significant" not in abstract.lower()


def test_manuscript_methods_lists_registered_algorithms_programmatically():
    methods = render_manuscript_methods(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    algorithm_names = list_algorithms()
    assert algorithm_names, "expected at least one registered algorithm"
    for name in algorithm_names:
        assert name in methods
    assert "Benjamini-Hochberg" in methods
    assert "Fisher's exact test" in methods
    assert "permutation test" in methods
    # n_genes_scanned reflects the actual length of the "genes" list passed in,
    # not the (larger, real-run) "genes_after_gating" count -- this fixture
    # only carries 2 sample gene rows.
    assert "2 scanned genes" in methods
    assert "q < 0.05" in methods


def test_manuscript_methods_uses_real_gating_and_config_counts():
    methods = render_manuscript_methods(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)
    assert "at least 5 distinct patient(s)" in methods
    assert "544 of 3919 genes passed the gate" in methods
    assert "4 using hand-curated gene configs" in methods
    assert "520 auto-configured" in methods


def test_manhattan_caption_counts_only_plottable_points_above_threshold():
    caption = render_manhattan_caption(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    # Both ETV6 and RET have a q-value and a composite score -> plottable.
    assert "2 scanned genes" in caption
    assert "q=0.05" in caption
    # Only ETV6 is FDR-significant among the plottable points.
    assert "1 gene above it" in caption


def test_manhattan_caption_handles_no_plottable_points():
    caption = render_manhattan_caption(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)
    assert "0 scanned genes" in caption
    assert "0 genes above it" in caption


def test_gene_highlight_significant_gene_matches_exact_text():
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][0]
    paragraph = render_gene_highlight(row)

    assert paragraph == (
        "ETV6 was analyzed across 90 fusion events, 71.1% in-frame and 75.6% domain-retained. "
        "Domain-retention Fisher's exact test p=5.12326e-06, genome-wide BH-adjusted "
        "q=0.00433427 (statistically significant at alpha=0.05)."
    )


def test_gene_highlight_honorable_mention_appends_note_verbatim():
    row = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["genes"][1]
    mention = PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS["honorable_mentions"][0]

    paragraph = render_gene_highlight(row, honorable_mention_note=mention["note"])

    assert paragraph.endswith(mention["note"])
    # RET's raw Fisher p-value (0.00042) is significant at alpha=0.05 even
    # though its FDR-adjusted q-value (0.12) is not -- the significance
    # clause is about the raw p-value, exactly like
    # cfh.reporting.text._domain_retention_paragraph already does.
    assert "statistically significant at alpha=0.05" in paragraph
    assert "not statistically significant" not in paragraph


def test_gene_highlight_states_missing_data_explicitly_not_invented():
    row = PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS["genes"][0]
    paragraph = render_gene_highlight(row)

    assert (
        "No domain-retention statistical test could be computed for FAKE1 in this scan."
        in paragraph
    )
    # Never invents a p-value/q-value that isn't there.
    assert "p=unavailable" not in paragraph


def test_discussion_bullets_are_gated_on_real_fields():
    bullets = render_discussion_bullets(PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS)

    joined = " ".join(bullets)
    assert "Benjamini-Hochberg FDR correction" in joined
    assert "conservative correction" in joined
    assert "2 scanned genes" in joined
    assert "4 gene(s) used hand-curated" in joined
    assert "520 gene(s) were auto-configured" in joined
    assert "kinase/catalytic-keyword heuristic" in joined
    assert "2026-09-04T14:36:53.295005+00:00" in joined
    assert "not validated clinical calls" in joined


def test_discussion_bullets_omit_auto_config_caveat_when_no_auto_configured_genes():
    payload = {**PAYLOAD_WITH_SIGNIFICANT_AND_HONORABLE_MENTIONS, "auto_config_gene_count": 0}
    bullets = render_discussion_bullets(payload)
    assert not any("auto-configured" in bullet for bullet in bullets)


def test_discussion_bullets_omit_generated_at_caveat_when_absent():
    bullets = render_discussion_bullets(PAYLOAD_WITH_NO_SIGNIFICANT_OR_HONORABLE_MENTIONS)
    assert not any("retrieved live" in bullet for bullet in bullets)
    # The standing clinical-use caveat is always present.
    assert any("not validated clinical calls" in bullet for bullet in bullets)

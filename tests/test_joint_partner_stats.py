import pytest

from cfh.model.fusion_event import FusionEvent
from cfh.stats.joint_partner_stats import (
    compute_pair_enrichment,
    evaluate_all_pairs,
    evaluate_gene_pair,
    extract_event_genes,
)


def test_extract_event_genes_directional():
    ev = FusionEvent(
        Event_id="e1",
        Cohort="c1",
        Five_prime_gene="EML4",
        Three_prime_gene="ALK",
    )
    assert extract_event_genes(ev, directional=True) == ("EML4", "ALK")
    assert extract_event_genes(ev, directional=False) == ("ALK", "EML4")


def test_extract_event_genes_fallback_to_site_genes():
    ev = FusionEvent(
        Event_id="e2",
        Cohort="c1",
        Site1_gene="EML4",
        Site2_gene="ALK",
    )
    assert extract_event_genes(ev, directional=True) == ("EML4", "ALK")
    assert extract_event_genes(ev, directional=False) == ("ALK", "EML4")


def test_extract_event_genes_missing_returns_none():
    ev = FusionEvent(Event_id="e3", Cohort="c1", Five_prime_gene="EML4")
    assert extract_event_genes(ev, directional=True) is None


def test_compute_pair_enrichment_fisher_overrepresented():
    # 100 events, 20 marginal 5p, 20 marginal 3p, expected = 4.0, observed = 18
    res = compute_pair_enrichment(
        observed_count=18,
        total_events=100,
        marginal_5p_count=20,
        marginal_3p_count=20,
        method="fisher",
        alpha=0.05,
    )
    assert res["p_value"] < 0.05
    assert res["is_significant"] is True
    assert res["expected_count"] == 4.0
    assert res["fold_enrichment"] == 4.5
    assert res["odds_ratio"] is not None and res["odds_ratio"] > 1.0


def test_compute_pair_enrichment_fisher_at_null():
    # 100 events, 20 marginal 5p, 20 marginal 3p, expected = 4.0, observed = 4
    res = compute_pair_enrichment(
        observed_count=4,
        total_events=100,
        marginal_5p_count=20,
        marginal_3p_count=20,
        method="fisher",
        alpha=0.05,
    )
    assert res["p_value"] >= 0.05
    assert res["is_significant"] is False
    assert res["expected_count"] == 4.0
    assert res["fold_enrichment"] == 1.0


def test_compute_pair_enrichment_binomial():
    res_sig = compute_pair_enrichment(
        observed_count=18,
        total_events=100,
        marginal_5p_count=20,
        marginal_3p_count=20,
        method="binomial",
        alpha=0.05,
    )
    assert res_sig["p_value"] < 0.05
    assert res_sig["is_significant"] is True

    res_null = compute_pair_enrichment(
        observed_count=4,
        total_events=100,
        marginal_5p_count=20,
        marginal_3p_count=20,
        method="binomial",
        alpha=0.05,
    )
    assert res_null["p_value"] >= 0.05
    assert res_null["is_significant"] is False


def test_compute_pair_enrichment_zero_events():
    res = compute_pair_enrichment(0, 0, 0, 0)
    assert res["p_value"] == 1.0
    assert res["is_significant"] is False
    assert res["expected_count"] == 0.0


def test_compute_pair_enrichment_invalid_method():
    with pytest.raises(ValueError, match="Unsupported method"):
        compute_pair_enrichment(5, 50, 10, 10, method="unknown_test")


def test_evaluate_gene_pair_and_all_pairs():
    events = [
        FusionEvent(
            Event_id=f"e{i}",
            Cohort="c1",
            Five_prime_gene="EML4",
            Three_prime_gene="ALK",
        )
        for i in range(10)
    ] + [
        FusionEvent(
            Event_id=f"o{i}",
            Cohort="c1",
            Five_prime_gene="OTHER_A",
            Three_prime_gene="OTHER_B",
        )
        for i in range(10)
    ]

    pair_res = evaluate_gene_pair(events, "EML4", "ALK")
    assert pair_res["gene_5p"] == "EML4"
    assert pair_res["gene_3p"] == "ALK"
    assert pair_res["observed_count"] == 10
    assert pair_res["total_events"] == 20
    assert pair_res["p_value"] < 0.05

    all_res = evaluate_all_pairs(events)
    assert len(all_res) == 2
    assert all_res[0]["gene_5p"] in ["EML4", "OTHER_A"]

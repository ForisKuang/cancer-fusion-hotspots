"""Tests for cohort-wide SV gene recurrence fetching and recurrence gating."""

from __future__ import annotations

from unittest.mock import MagicMock

from cfh.cohort.recurrence import (
    DEFAULT_MIN_DISTINCT_PATIENTS,
    GeneRecurrence,
    fetch_cohort_gene_recurrence,
    gate_genes_by_recurrence,
)
from cfh.ingestion import cbioportal_api


def test_fetch_structural_variant_genes_posts_expected_body_without_real_network():
    mock_session = MagicMock()
    mock_response = MagicMock(status_code=200)
    mock_response.json.return_value = [
        {
            "hugoGeneSymbol": "BRAF",
            "entrezGeneId": 673,
            "numberOfAlteredCases": 247,
            "totalCount": 251,
        }
    ]
    mock_session.post.return_value = mock_response

    result = cbioportal_api.fetch_structural_variant_genes(
        ["msk_impact_50k_2026"], session=mock_session
    )

    assert result == mock_response.json.return_value
    mock_session.post.assert_called_once()
    called_url = mock_session.post.call_args.args[0]
    _, kwargs = mock_session.post.call_args
    assert called_url == f"{cbioportal_api.DEFAULT_BASE_URL}/structuralvariant-genes/fetch"
    assert kwargs["json"] == {"studyIds": ["msk_impact_50k_2026"]}


def test_fetch_structural_variant_genes_retries_transient_failure():
    mock_session = MagicMock()
    unavailable = MagicMock(status_code=503)
    recovered = MagicMock(status_code=200)
    recovered.json.return_value = []
    mock_session.post.side_effect = [unavailable, recovered]

    result = cbioportal_api.fetch_structural_variant_genes(
        ["some_study"], session=mock_session, max_retries=1, backoff_seconds=0
    )

    assert result == []
    assert mock_session.post.call_count == 2


def test_fetch_cohort_gene_recurrence_parses_alteration_counts():
    mock_session = MagicMock()
    mock_response = MagicMock(status_code=200)
    mock_response.json.return_value = [
        {
            "hugoGeneSymbol": "BRAF",
            "entrezGeneId": 673,
            "numberOfAlteredCases": 247,
            "totalCount": 251,
        },
        {
            "hugoGeneSymbol": "SINGLETON1",
            "entrezGeneId": 111,
            "numberOfAlteredCases": 1,
            "totalCount": 1,
        },
        # A malformed record with no gene symbol must be skipped, not crash the parse.
        {"entrezGeneId": 999, "numberOfAlteredCases": 3, "totalCount": 3},
    ]
    mock_session.post.return_value = mock_response

    recurrence = fetch_cohort_gene_recurrence("msk_impact_50k_2026", session=mock_session)

    assert len(recurrence) == 2
    braf = next(gene for gene in recurrence if gene.hugo_gene_symbol == "BRAF")
    assert braf.entrez_gene_id == 673
    assert braf.distinct_patient_count == 247
    assert braf.total_sv_count == 251


def test_gate_genes_by_recurrence_reports_full_count_before_filtering():
    recurrence = [
        GeneRecurrence("BRAF", 673, 247, 251),
        GeneRecurrence("RET", 5979, 219, 230),
        GeneRecurrence("SINGLETON1", 1, 1, 1),
        GeneRecurrence("SINGLETON2", 2, 1, 1),
        GeneRecurrence("BORDERLINE", 3, 5, 5),
    ]

    gate = gate_genes_by_recurrence(recurrence, min_distinct_patients=5, study_id="some_study")

    # The pre-filter universe must always be visible, not just the passing subset.
    assert gate.total_genes == 5
    assert gate.passing_count == 3
    assert gate.filtered_out_count == 2
    assert {gene.hugo_gene_symbol for gene in gate.passing_genes} == {"BRAF", "RET", "BORDERLINE"}
    assert {gene.hugo_gene_symbol for gene in gate.filtered_out_genes} == {
        "SINGLETON1",
        "SINGLETON2",
    }
    # Every input gene ends up in exactly one bucket.
    assert gate.passing_count + gate.filtered_out_count == gate.total_genes


def test_gate_genes_by_recurrence_default_threshold_is_five():
    assert DEFAULT_MIN_DISTINCT_PATIENTS == 5
    recurrence = [GeneRecurrence("A", 1, 4, 4), GeneRecurrence("B", 2, 5, 5)]
    gate = gate_genes_by_recurrence(recurrence)
    assert [gene.hugo_gene_symbol for gene in gate.passing_genes] == ["B"]


def test_gate_genes_by_recurrence_sorts_passing_genes_by_descending_recurrence():
    recurrence = [
        GeneRecurrence("LOW", 1, 6, 6),
        GeneRecurrence("HIGH", 2, 500, 500),
        GeneRecurrence("MID", 3, 50, 50),
    ]
    gate = gate_genes_by_recurrence(recurrence, min_distinct_patients=5)
    assert [gene.hugo_gene_symbol for gene in gate.passing_genes] == ["HIGH", "MID", "LOW"]


def test_gate_genes_by_recurrence_rejects_negative_threshold():
    import pytest

    with pytest.raises(ValueError):
        gate_genes_by_recurrence([], min_distinct_patients=-1)

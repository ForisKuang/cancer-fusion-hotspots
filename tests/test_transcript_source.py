from unittest.mock import MagicMock

import pytest

from cfh.genes.registry import load_gene_config
from cfh.mapping import transcript_source


def test_annotation_with_explicit_exon_number_is_marked_source_annotation():
    gene_config = load_gene_config("braf")
    result = transcript_source.resolve_transcript_mapping(
        "KIAA1549-BRAF fusion, exon16:exon9, in-frame", gene_config
    )
    assert result.source == "annotation"
    assert result.breakpoint_exon == 16


def test_annotation_lacking_exon_info_falls_back_to_ensembl(monkeypatch):
    gene_config = load_gene_config("braf")
    mock_client = MagicMock()
    mock_client.fetch_protein_features.return_value = [{"type": "domain"}]

    result = transcript_source.resolve_transcript_mapping(
        None,
        gene_config,
        ensembl_protein_id="ENSP00000288602",
        ensembl_client=mock_client,
    )

    assert result.source == "ensembl_fallback"
    assert result.breakpoint_exon is None
    mock_client.fetch_protein_features.assert_called_once_with("ENSP00000288602")


def test_fallback_with_breakpoint_coordinate_returns_real_overlapping_features():
    """The fallback must actually use the breakpoint, not just fetch and ignore it."""
    gene_config = load_gene_config("braf")
    mock_client = MagicMock()
    mock_client.fetch_protein_features.return_value = [
        {"type": "protein_feature", "id": "PF07714", "description": "Protein kinase domain",
         "start": 457, "end": 717},
        {"type": "protein_feature", "id": "PF00130", "description": "RAS-binding domain",
         "start": 155, "end": 227},
    ]

    result = transcript_source.resolve_transcript_mapping(
        None,
        gene_config,
        breakpoint_aa=500,
        ensembl_protein_id="ENSP00000288602",
        ensembl_client=mock_client,
    )

    assert result.source == "ensembl_fallback"
    assert result.breakpoint_exon is None  # genuinely unavailable from this data source
    assert result.breakpoint_protein_features is not None
    assert len(result.breakpoint_protein_features) == 1
    assert result.breakpoint_protein_features[0]["id"] == "PF07714"


def test_fallback_unavailable_when_no_ensembl_protein_id_supplied():
    gene_config = load_gene_config("braf")
    with pytest.raises(transcript_source.EnsemblFallbackUnavailable):
        transcript_source.resolve_transcript_mapping(
            None,
            gene_config,
            breakpoint_aa=500,
        )


def test_ensembl_client_calls_expected_url_and_params(monkeypatch):
    mock_session = MagicMock()
    mock_session.get.return_value.json.return_value = []
    client = transcript_source.EnsemblClient(session=mock_session)

    client.fetch_protein_features("ENSP00000288602")

    args, kwargs = mock_session.get.call_args
    assert args[0] == "https://rest.ensembl.org/overlap/translation/ENSP00000288602"
    assert kwargs["params"] == {"feature": "protein_feature"}

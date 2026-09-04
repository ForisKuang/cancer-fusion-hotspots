import json
from unittest.mock import MagicMock

import pytest

from cfh.genes.registry import load_gene_config
from cfh.mapping import genome_nexus_source as gns
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
        {
            "type": "protein_feature",
            "id": "PF07714",
            "description": "Protein kinase domain",
            "start": 457,
            "end": 717,
        },
        {
            "type": "protein_feature",
            "id": "PF00130",
            "description": "RAS-binding domain",
            "start": 155,
            "end": 227,
        },
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


def test_resolve_breakpoint_protein_position_prefers_annotation_exon():
    gene_config = load_gene_config("braf")
    result = transcript_source.resolve_breakpoint_protein_position(
        "KIAA1549-BRAF fusion, exon16:exon9, in-frame", gene_config
    )
    assert result.source == "annotation"
    assert result.breakpoint_exon == 16


def test_resolve_breakpoint_protein_position_uses_genome_nexus_for_real_estimate(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Regression for the CDS/UTR-bounds bug: the PRODUCTION entry point
    (not the lower-level CDS-offset helper in isolation) must derive CDS
    bounds from the transcript's UTRs and use them, or it silently counts
    5' UTR nucleotides as coding sequence. For this exact fixture and
    breakpoint, the correct residue is 217; the pre-fix bug (ignoring UTR
    bounds, counting 61nt of 5' UTR as CDS) returned 237.
    """
    gene_config = load_gene_config("braf")
    payload = json.loads(genome_nexus_canonical_transcript_fixture_path.read_text())
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.return_value = payload

    canonical = gns.parse_canonical_transcript(payload)
    exon5 = next(e for e in canonical.exons if e.rank == 5)
    breakpoint_genomic = exon5.start + 62

    result = transcript_source.resolve_breakpoint_protein_position(
        None,
        gene_config,
        breakpoint_genomic=breakpoint_genomic,
        genome_nexus_client=mock_gn_client,
    )

    assert result.source == "genome_nexus_fallback"
    assert result.is_intronic_breakpoint is False
    assert result.breakpoint_exon == 5
    assert result.transcript_id == "NM_004333"
    assert result.breakpoint_protein_position == 217
    assert result.breakpoint_protein_position != 237  # the pre-fix (UTR-ignoring) value


def test_resolve_breakpoint_protein_position_falls_back_to_ensembl_with_protein_coordinate():
    """When Genome Nexus doesn't know the gene but the caller already has a
    protein-level breakpoint_aa, the Ensembl protein-feature fallback can
    still produce a real (non-empty) result.
    """
    gene_config = load_gene_config("braf")
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.side_effect = gns.GenomeNexusGeneNotFound("nope")
    mock_ensembl_client = MagicMock()
    mock_ensembl_client.fetch_protein_features.return_value = [
        {"id": "PF07714", "start": 457, "end": 717},
    ]

    result = transcript_source.resolve_breakpoint_protein_position(
        None,
        gene_config,
        breakpoint_genomic=140507800,
        breakpoint_aa=500,
        genome_nexus_client=mock_gn_client,
        ensembl_protein_id="ENSP00000288602",
        ensembl_client=mock_ensembl_client,
    )

    assert result.source == "ensembl_fallback"
    assert result.breakpoint_protein_features is not None
    assert len(result.breakpoint_protein_features) == 1
    mock_ensembl_client.fetch_protein_features.assert_called_once_with("ENSP00000288602")


def test_resolve_breakpoint_protein_position_uses_canonical_protein_id_without_redundant_arg(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Regression: when Genome Nexus finds the gene (so its canonical
    payload's protein_id is already known) and a real breakpoint_aa is
    supplied, the Ensembl fallback must use that protein_id automatically
    -- the caller should not have to redundantly pass ensembl_protein_id
    just to unlock a mapping the code already has enough information for.
    """
    gene_config = load_gene_config("braf")
    payload = json.loads(genome_nexus_canonical_transcript_fixture_path.read_text())
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.return_value = payload

    mock_ensembl_client = MagicMock()
    mock_ensembl_client.fetch_protein_features.return_value = [
        {"id": "PF07714", "start": 457, "end": 717},
    ]

    result = transcript_source.resolve_breakpoint_protein_position(
        None,
        gene_config,
        breakpoint_aa=500,
        genome_nexus_client=mock_gn_client,
        ensembl_client=mock_ensembl_client,
        # Deliberately no ensembl_protein_id -- it must come from the
        # canonical transcript Genome Nexus already returned.
    )

    assert result.source == "ensembl_fallback"
    assert result.breakpoint_protein_features is not None
    assert len(result.breakpoint_protein_features) == 1
    mock_ensembl_client.fetch_protein_features.assert_called_once_with("ENSP00000288602")


def test_resolve_breakpoint_protein_position_raises_rather_than_silent_empty_success():
    """A genomic-only breakpoint with no Genome Nexus mapping must raise --
    NOT return an apparently-successful ensembl_fallback with no computed
    position, since the Ensembl protein-feature endpoint cannot consume a
    genomic coordinate at all.
    """
    gene_config = load_gene_config("braf")
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.side_effect = gns.GenomeNexusGeneNotFound("nope")
    mock_ensembl_client = MagicMock()

    with pytest.raises(transcript_source.EnsemblFallbackUnavailable):
        transcript_source.resolve_breakpoint_protein_position(
            None,
            gene_config,
            breakpoint_genomic=140507800,
            genome_nexus_client=mock_gn_client,
            ensembl_protein_id="ENSP00000288602",
            ensembl_client=mock_ensembl_client,
        )

    mock_ensembl_client.fetch_protein_features.assert_not_called()


def test_resolve_breakpoint_protein_position_raises_when_genuinely_unmappable():
    gene_config = load_gene_config("braf")
    mock_gn_client = MagicMock()
    mock_gn_client.fetch_canonical_transcript.side_effect = gns.GenomeNexusGeneNotFound("nope")

    with pytest.raises(transcript_source.EnsemblFallbackUnavailable):
        transcript_source.resolve_breakpoint_protein_position(
            None,
            gene_config,
            breakpoint_genomic=140507800,
            genome_nexus_client=mock_gn_client,
        )

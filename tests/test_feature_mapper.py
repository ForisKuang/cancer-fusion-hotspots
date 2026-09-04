import json
from unittest.mock import MagicMock

import pytest

from cfh.genes.registry import GeneConfig, KeyDomain, load_gene_config
from cfh.mapping import feature_mapper
from cfh.mapping.domain_source import ProteinDomain, UniProtDomainSource
from cfh.model.fusion_event import FusionEvent


@pytest.fixture(autouse=True)
def _reset_default_domain_source():
    feature_mapper.reset_default_domain_source()
    yield
    feature_mapper.reset_default_domain_source()


def _domain_source(uniprot_fixture_path):
    payload = json.loads(uniprot_fixture_path.read_text())
    source = UniProtDomainSource()
    source._cache["P15056"] = source.parse(payload)
    return source


def _kinase_bounds(uniprot_fixture_path):
    payload = json.loads(uniprot_fixture_path.read_text())
    feature = next(
        f
        for f in payload["features"]
        if f.get("type") == "Domain" and "kinase" in (f.get("description") or "").lower()
    )
    return feature["location"]["start"]["value"], feature["location"]["end"]["value"]


def _event(sample_id: str) -> FusionEvent:
    return FusionEvent(
        Event_id=f"EVT-{sample_id}", Cohort="msk_impact_50k_2026", Sample_id=sample_id
    )


def test_breakpoint_after_full_kinase_domain_is_retained(uniprot_fixture_path):
    gene_config = load_gene_config("braf")
    domain_source = _domain_source(uniprot_fixture_path)
    _, kinase_end = _kinase_bounds(uniprot_fixture_path)

    feature = feature_mapper.map_event(
        _event("SAMPLE-A"),
        gene_config,
        role="five_prime",
        junction_position_aa=kinase_end + 10,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags["kinase"] == "retained"


def test_breakpoint_inside_kinase_domain_is_disrupted(uniprot_fixture_path):
    gene_config = load_gene_config("braf")
    domain_source = _domain_source(uniprot_fixture_path)
    kinase_start, kinase_end = _kinase_bounds(uniprot_fixture_path)
    midpoint = (kinase_start + kinase_end) // 2

    feature = feature_mapper.map_event(
        _event("SAMPLE-B"),
        gene_config,
        role="five_prime",
        junction_position_aa=midpoint,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags["kinase"] == "disrupted"


def test_breakpoint_before_kinase_domain_is_lost(uniprot_fixture_path):
    gene_config = load_gene_config("braf")
    domain_source = _domain_source(uniprot_fixture_path)
    kinase_start, _ = _kinase_bounds(uniprot_fixture_path)

    feature = feature_mapper.map_event(
        _event("SAMPLE-C"),
        gene_config,
        role="five_prime",
        junction_position_aa=kinase_start - 50,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags["kinase"] == "lost"


def test_breakpoint_exactly_at_domain_end_is_retained(uniprot_fixture_path):
    gene_config = load_gene_config("braf")
    domain_source = _domain_source(uniprot_fixture_path)
    _, kinase_end = _kinase_bounds(uniprot_fixture_path)

    feature = feature_mapper.map_event(
        _event("SAMPLE-D"),
        gene_config,
        role="five_prime",
        junction_position_aa=kinase_end,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags["kinase"] == "retained"


def test_breakpoint_exactly_at_domain_start_is_disrupted(uniprot_fixture_path):
    gene_config = load_gene_config("braf")
    domain_source = _domain_source(uniprot_fixture_path)
    kinase_start, _ = _kinase_bounds(uniprot_fixture_path)

    feature = feature_mapper.map_event(
        _event("SAMPLE-E"),
        gene_config,
        role="five_prime",
        junction_position_aa=kinase_start,
        domain_source=domain_source,
    )

    # Only the single boundary residue would be retained -- the rest of the
    # domain is cut, so this is a disruption, not a clean retain/lose split.
    assert feature.Domain_retention_flags["kinase"] == "disrupted"


def test_disruption_required_domains_get_retention_flags_alongside_key_domains():
    """BRAF's opt-in disruption_required_domains (RAS-binding, cysteine-rich)
    are populated through the same map_event pass as key_domains, using real
    Genome Nexus Pfam accessions/boundaries -- no separate mapping path."""
    gene_config = load_gene_config("braf")
    domain_source = MagicMock()
    domain_source.fetch.return_value = [
        ProteinDomain(
            name="PF07714", start_aa=458, end_aa=712, source="genome_nexus", accession="PF07714"
        ),
        ProteinDomain(
            name="PF02196", start_aa=156, end_aa=227, source="genome_nexus", accession="PF02196"
        ),
        ProteinDomain(
            name="PF00130", start_aa=235, end_aa=280, source="genome_nexus", accession="PF00130"
        ),
    ]

    # Breakpoint at aa 380, BRAF contributing its 3' (C-terminal) fragment:
    # the kinase domain (458-712) is downstream and retained, while both
    # N-terminal regulatory domains (156-227, 235-280) are excised.
    feature = feature_mapper.map_event(
        _event("SAMPLE-REAL"),
        gene_config,
        role="three_prime",
        junction_position_aa=380,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags["kinase"] == "retained"
    assert feature.Domain_retention_flags["ras_binding"] == "lost"
    assert feature.Domain_retention_flags["cysteine_rich"] == "lost"
    assert "Protein kinase domain" in feature.Retained_domains
    assert "RAS-binding domain" in feature.Lost_domains
    assert "Cysteine-rich domain" in feature.Lost_domains


def test_domain_listed_in_both_key_and_disruption_domains_is_not_double_processed():
    """Regression: a gene that (however unusually) lists the same domain
    key in both key_domains and disruption_required_domains must have it
    mapped/classified exactly once, not twice -- key_domains' entry wins,
    since it's listed first."""
    gene_config = GeneConfig(
        gene_symbol="FAKE3",
        canonical_transcript_id="NM_000003",
        protein_id="P00003",
        key_domains=[KeyDomain(name="Shared domain (retention copy)", source="test", key="shared")],
        disruption_required_domains=[
            KeyDomain(name="Shared domain (disruption copy)", source="test", key="shared")
        ],
    )
    domain_source = MagicMock()
    domain_source.fetch.return_value = [
        ProteinDomain(name="shared", start_aa=100, end_aa=200, source="test", accession=None)
    ]

    # Breakpoint downstream of the domain, five_prime role: the domain
    # (100-200) falls fully inside the retained N-terminal fragment.
    feature = feature_mapper.map_event(
        _event("SAMPLE-DUP"),
        gene_config,
        role="five_prime",
        junction_position_aa=300,
        domain_source=domain_source,
    )

    assert feature.Domain_retention_flags == {"shared": "retained"}
    assert feature.Retained_domains == ["Shared domain (retention copy)"]
    assert feature.Lost_domains == []
    assert feature.Disrupted_domains == []


def test_default_domain_source_is_shared_across_calls_without_explicit_source(
    uniprot_fixture_path,
):
    """The real call path (no manually-shared domain_source) must still cache."""
    gene_config = load_gene_config("braf")
    mock_session = MagicMock()
    mock_session.get.return_value.json.return_value = json.loads(uniprot_fixture_path.read_text())
    mock_session.get.return_value.raise_for_status.return_value = None

    shared_source = UniProtDomainSource(session=mock_session)
    feature_mapper._default_domain_source = shared_source

    feature_mapper.map_event(
        _event("SAMPLE-F"), gene_config, role="five_prime", junction_position_aa=600
    )
    feature_mapper.map_event(
        _event("SAMPLE-G"), gene_config, role="five_prime", junction_position_aa=400
    )

    mock_session.get.assert_called_once()

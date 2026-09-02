import json
from unittest.mock import MagicMock

import pytest

from cfh.genes.registry import load_gene_config
from cfh.mapping import feature_mapper
from cfh.mapping.domain_source import UniProtDomainSource
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


def test_default_domain_source_is_shared_across_calls_without_explicit_source(
    uniprot_fixture_path,
):
    """The real call path (no manually-shared domain_source) must still cache."""
    gene_config = load_gene_config("braf")
    mock_session = MagicMock()
    mock_session.get.return_value.json.return_value = json.loads(
        uniprot_fixture_path.read_text()
    )
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

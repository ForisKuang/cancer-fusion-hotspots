"""CI-gating BRAF kinase-domain-retention benchmark for MSK-IMPACT 50k."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.feature_mapper import map_event
from cfh.mapping.transcript_source import resolve_breakpoint_protein_position
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import fishers_frame_domain_test, permutation_null_test

_FIXTURE = Path(__file__).parents[1] / "fixtures" / "braf_fusion_calls_benchmark.json"
_PROFILE_ID = "msk_impact_50k_2026_structural_variants"


def _fixture_events_and_features() -> tuple[list[FusionEvent], list[FusionFeature]]:
    # This fixture is synthetic, modeled on the named BRAF-fusion examples in
    # Zehir et al. (PMC5461196), pending a real MSK-IMPACT 50k pull through
    # the existing cBioPortal client.  It deliberately includes non-recurrent
    # out-of-frame/kinase-lost comparison events for the formal test.
    payload = json.loads(_FIXTURE.read_text())
    events = []
    features = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort="msk_impact_50k_2026",
                Sample_id=row["event_id"],
                Fusion_name=row["fusion_name"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=row["is_protein_fusion"],
                Three_prime_gene="BRAF",
            )
        )
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene="BRAF",
                Role="three_prime",
                Transcript_id="NM_004333",
                Junction_position_aa=row["junction_position_aa"],
                Domain_retention_flags={"kinase": row["kinase_status"]},
            )
        )
    return events, features


def _assert_benchmark(events: list[FusionEvent], features: list[FusionFeature]) -> None:
    config = load_gene_config("braf")
    result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 1_000}
    )
    assert result.Summary["fisher_p_value"] < 0.05


def test_braf_kinase_retention_is_enriched_in_fixture():
    events, features = _fixture_events_and_features()
    _assert_benchmark(events, features)


def test_fisher_handles_the_paper_like_zero_off_diagonal_cell():
    odds_ratio, p_value = fishers_frame_domain_test([[6, 0], [0, 8]])
    assert odds_ratio == float("inf")
    assert 0 <= p_value <= 1


def test_permutation_null_is_bit_identical_for_a_fixed_seed():
    events, features = _fixture_events_and_features()
    config = load_gene_config("braf")
    first = permutation_null_test(events, features, config, seed=1729, n_permutations=250)
    second = permutation_null_test(events, features, config, seed=1729, n_permutations=250)
    assert first == second


def _call_value(call: dict, snake_name: str, camel_name: str):
    return call.get(snake_name, call.get(camel_name))


def _real_events_and_features(calls: list[dict]) -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Map real cBioPortal calls through the reviewed breakpoint/domain path."""
    config = load_gene_config("braf")
    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for index, call in enumerate(calls):
        site1 = _call_value(call, "site1_hugo_symbol", "site1HugoSymbol")
        site2 = _call_value(call, "site2_hugo_symbol", "site2HugoSymbol")
        connection = _call_value(call, "connection_type", "connectionType")
        if site1 != "BRAF" and site2 != "BRAF":
            continue
        if connection == "5to3" and site2 == "BRAF":
            role, breakpoint = "three_prime", _call_value(call, "site2_position", "site2Position")
        elif connection == "3to5" and site1 == "BRAF":
            role, breakpoint = "three_prime", _call_value(call, "site1_position", "site1Position")
        elif connection == "5to3" and site1 == "BRAF":
            role, breakpoint = "five_prime", _call_value(call, "site1_position", "site1Position")
        elif connection == "3to5" and site2 == "BRAF":
            role, breakpoint = "five_prime", _call_value(call, "site2_position", "site2Position")
        else:
            continue
        if breakpoint is None:
            continue
        frame = _call_value(call, "site2_effect_on_frame", "site2EffectOnFrame")
        event = FusionEvent(
            Event_id=str(
                _call_value(call, "structural_variant_id", "structuralVariantId") or index
            ),
            Cohort="msk_impact_50k_2026",
            Site1_gene=site1,
            Site2_gene=site2,
            Frame_status=frame,
            Is_protein_fusion=True,
        )
        try:
            mapping = resolve_breakpoint_protein_position(
                None, config, breakpoint_genomic=int(breakpoint)
            )
        except Exception:
            continue
        events.append(event)
        features.append(
            map_event(
                event,
                config,
                role=role,
                junction_position_aa=mapping.breakpoint_protein_position,
            )
        )
    return events, features


@pytest.mark.network
@pytest.mark.skipif(
    os.getenv("CFH_RUN_NETWORK_TESTS") != "1",
    reason="set CFH_RUN_NETWORK_TESTS=1 to run live cBioPortal/Genome Nexus benchmark",
)
def test_braf_kinase_retention_in_real_msk_impact_50k():
    calls = cbioportal_api.fetch_structural_variants(
        entrez_gene_ids=[673], molecular_profile_ids=[_PROFILE_ID]
    )
    events, features = _real_events_and_features(calls)
    _assert_benchmark(events, features)

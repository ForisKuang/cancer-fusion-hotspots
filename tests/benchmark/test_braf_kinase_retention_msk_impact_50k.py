"""CI-gating BRAF kinase-domain-retention benchmark for MSK-IMPACT 50k."""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from cfh.algorithms.domain_retention import DomainRetentionAlgorithm
from cfh.algorithms.frequency import FrequencyAnalysis
from cfh.genes.registry import load_gene_config
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.real_benchmark import analyze_structural_variant_calls
from cfh.stats.breakpoint_tests import (
    build_frame_domain_contingency_table,
    fishers_frame_domain_test,
    permutation_null_test,
)

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
    config = load_gene_config("braf")
    table = build_frame_domain_contingency_table(events, features, config)
    odds_ratio, p_value = fishers_frame_domain_test(table)

    assert table == [[6, 0], [0, 8]]
    assert odds_ratio == float("inf")
    assert p_value == pytest.approx(0.000333000333000333)

    result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 1_000}
    )
    assert result.Tables["frame_domain_contingency_table"] == table
    assert result.Summary["fisher_odds_ratio"] == odds_ratio
    assert result.Summary["fisher_p_value"] == pytest.approx(p_value)


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


def _real_events_and_features(
    calls: list[dict], genome_nexus_client: GenomeNexusClient | None = None
) -> tuple[list[FusionEvent], list[FusionFeature]]:
    """Use the production API adapter, normalizer, and Genome Nexus mapping path."""
    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=genome_nexus_client,
        n_permutations=25,
    )
    return run.events, run.features


def test_real_call_field_mapping_uses_production_normalizer_and_genome_nexus(
    genome_nexus_canonical_transcript_fixture_path,
):
    """Regression: raw ``NA`` frame values must not bypass production normalization."""
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(
        genome_nexus_canonical_transcript_fixture_path.read_text()
    )
    calls = [
        {
            "sampleId": "SAMPLE-1",
            "site1HugoSymbol": "PARTNER1",
            "site2HugoSymbol": "BRAF",
            "connectionType": "3to3",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "eventInfo": "Protein Fusion: in frame  {PARTNER1:BRAF}",
        },
        {
            "sampleId": "SAMPLE-2",
            "site1HugoSymbol": "BRAF",
            "site2HugoSymbol": "PARTNER2",
            "connectionType": "5to5",
            "site1Position": 140493152,
            "site2EffectOnFrame": "NA",
            "eventInfo": "Protein Fusion: out of frame  {BRAF:PARTNER2}",
        },
        {
            "sampleId": "SAMPLE-3",
            "site1HugoSymbol": "BRAF",
            "site2HugoSymbol": "PARTNER3",
            "connectionType": "5to3",
            "site1Position": 140493152,
            "site2EffectOnFrame": "NA",
            "eventInfo": "Protein Fusion: mid-exon  {PARTNER3:BRAF}",
        },
        {
            "sampleId": "SAMPLE-4",
            "site1HugoSymbol": "PARTNER4",
            "site2HugoSymbol": "BRAF",
            "connectionType": "3to5",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "eventInfo": "Protein Fusion: in frame  {PARTNER4:BRAF}",
        },
    ]

    events, features = _real_events_and_features(calls, client)

    assert [event.Frame_status for event in events] == [
        "in-frame",
        "out-of-frame",
        "unknown",
        "in-frame",
    ]
    assert [feature.Role for feature in features] == [
        "three_prime",
        "five_prime",
        "three_prime",
        "three_prime",
    ]
    assert all(feature.Junction_position_aa is not None for feature in features)
    assert all(feature.Domain_retention_flags["kinase"] != "unknown" for feature in features)
    client.fetch_canonical_transcript.assert_called()


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
    config = load_gene_config("braf")
    result = DomainRetentionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 1_000}
    )
    frequency = FrequencyAnalysis().run(events, features, config, {})

    assert len(calls) > 0
    assert len(events) == len(features) > 0
    assert any(event.Frame_status == "in-frame" for event in events)
    assert all(feature.Domain_retention_flags["kinase"] != "unknown" for feature in features)
    assert 0 <= result.Summary["fisher_p_value"] <= 1
    assert sum(row["Event_count"] for row in frequency.Tables["Partner_gene_counts"]) == len(
        events
    )
    print(
        json.dumps(
            {
                "raw_structural_variants": len(calls),
                "protein_fusions": len(events),
                "in_frame": sum(event.Frame_status == "in-frame" for event in events),
                "kinase_retained": sum(
                    feature.Domain_retention_flags["kinase"] == "retained"
                    for feature in features
                ),
                "fisher_p_value": result.Summary["fisher_p_value"],
            },
            sort_keys=True,
        )
    )

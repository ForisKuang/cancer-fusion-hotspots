from cfh.algorithms import registry
from cfh.algorithms.cutpoint_detection import CutpointDetectionAlgorithm
from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

# A fabricated single-gene config, deliberately unrelated to any real gene
# in genes/configs/, to demonstrate the algorithm has no gene-specific
# hardcoding of its own.
_FAKE_GENE = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[KeyDomain(name="Made-up domain", source="test", key="made_up")],
)


def _events_and_features(records: list[tuple[str, int, str]]) -> tuple[list, list]:
    events = []
    features = []
    for event_id, position, status in records:
        events.append(FusionEvent(Event_id=event_id, Cohort="synthetic"))
        features.append(
            FusionFeature(
                Event_id=event_id,
                Gene="FAKE1",
                Junction_position_aa=position,
                Domain_retention_flags={"made_up": status},
            )
        )
    return events, features


def _clean_separation_records() -> list[tuple[str, int, str]]:
    retained = [("r0", 100, "retained"), ("r1", 110, "retained"), ("r2", 120, "retained")]
    lost = [("l0", 300, "lost"), ("l1", 310, "disrupted"), ("l2", 320, "lost")]
    return retained + lost


def test_algorithm_registered():
    assert registry.get("cutpoint_detection") is CutpointDetectionAlgorithm
    assert "cutpoint_detection" in registry.list_algorithms()


def test_recovers_cutpoint_from_synthetic_clean_separation():
    events, features = _events_and_features(_clean_separation_records())

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 300}
    )

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "cutpoint_detection"
    assert result.Summary["determinable"] is True
    assert result.Summary["inferred_cutpoint_aa"] == 120
    assert result.Summary["n_events_analyzed"] == 6
    assert 0.0 <= result.Summary["corrected_p_value"] <= 1.0
    assert result.Warnings == []
    assert result.Tables["cutpoint_scan"]
    assert result.Summary["known_domain_boundary_comparison"] is None


def test_domain_boundary_comparison_reports_nearest_distance():
    events, features = _events_and_features(_clean_separation_records())

    result = CutpointDetectionAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {"seed": 42, "n_permutations": 300, "domain_boundaries": [118, 400]},
    )

    comparison = result.Summary["known_domain_boundary_comparison"]
    assert comparison == {"nearest_known_domain_boundary_aa": 118, "distance_aa": 2}


def test_degenerate_input_is_reported_as_not_determinable_without_raising():
    events, features = _events_and_features(
        [("a", 100, "retained"), ("b", 110, "retained"), ("c", 120, "retained")]
    )

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["determinable"] is False
    assert result.Summary["inferred_cutpoint_aa"] is None
    assert result.Warnings
    assert result.Summary["known_domain_boundary_comparison"] is None


def test_events_missing_mapped_positions_are_excluded_not_crashed_on():
    events = [FusionEvent(Event_id="unmapped", Cohort="synthetic")]
    features = [
        FusionFeature(
            Event_id="unmapped",
            Gene="FAKE1",
            Junction_position_aa=None,
            Domain_retention_flags={"made_up": "retained"},
        )
    ]

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["n_events_analyzed"] == 0
    assert result.Summary["determinable"] is False


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _events_and_features(_clean_separation_records())
    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)

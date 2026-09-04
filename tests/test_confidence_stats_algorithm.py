import pytest

from cfh.algorithms import registry
from cfh.algorithms.confidence_stats import ConfidenceStatsAlgorithm
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent


def _make_events() -> list[FusionEvent]:
    events = []
    # Group A: Is_protein_fusion=True, higher read support, mostly antisense.
    for i in range(6):
        events.append(
            FusionEvent(
                Event_id=f"a{i}",
                Cohort="cohort1",
                Is_protein_fusion=True,
                Is_antisense=(i < 5),
                Paired_end_read_support=50 + i,
            )
        )
    # Group B: Is_protein_fusion=False, lower read support, rarely antisense.
    for i in range(6):
        events.append(
            FusionEvent(
                Event_id=f"b{i}",
                Cohort="cohort1",
                Is_protein_fusion=False,
                Is_antisense=(i < 1),
                Paired_end_read_support=5 + i,
            )
        )
    return events


def test_algorithm_registered():
    assert registry.get("confidence_stats") is ConfidenceStatsAlgorithm
    assert "confidence_stats" in registry.list_algorithms()


def test_both_mle_and_ttest_requested():
    events = _make_events()
    params = {
        "group_field": "Is_protein_fusion",
        "outcome_field": "Is_antisense",
        "numeric_field": "Paired_end_read_support",
    }

    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "confidence_stats"
    assert result.Summary is not None
    assert "mle" in result.Summary
    assert "ttest" in result.Summary

    mle = result.Summary["mle"]
    assert mle["outcome_field"] == "Is_antisense"
    assert set(mle["groups"].keys()) == {"True", "False"}
    for group_stats in mle["groups"].values():
        low = group_stats["ci_low"]
        estimate = group_stats["point_estimate"]
        high = group_stats["ci_high"]
        assert 0.0 <= low <= estimate <= high <= 1.0

    ttest = result.Summary["ttest"]
    assert ttest["numeric_field"] == "Paired_end_read_support"
    assert ttest["p_value"] < 0.05  # groups are clearly separated by construction


def test_mle_only_requested_omits_ttest_block():
    events = _make_events()
    params = {"group_field": "Is_protein_fusion", "outcome_field": "Is_antisense"}

    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert "mle" in result.Summary
    assert result.Summary.get("ttest") is None
    assert "ttest" not in result.Summary


def test_ttest_only_requested_omits_mle_block():
    events = _make_events()
    params = {"group_field": "Is_protein_fusion", "numeric_field": "Paired_end_read_support"}

    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert "ttest" in result.Summary
    assert result.Summary.get("mle") is None
    assert "mle" not in result.Summary


def test_requires_group_field():
    events = _make_events()
    with pytest.raises(ValueError):
        ConfidenceStatsAlgorithm().run(events, [], None, {"outcome_field": "Is_antisense"})


def test_requires_at_least_one_of_outcome_or_numeric_field():
    events = _make_events()
    with pytest.raises(ValueError):
        ConfidenceStatsAlgorithm().run(events, [], None, {"group_field": "Is_protein_fusion"})


def test_result_schema_matches_canonical_algorithm_result_fields():
    events = _make_events()
    params = {
        "group_field": "Is_protein_fusion",
        "outcome_field": "Is_antisense",
        "numeric_field": "Paired_end_read_support",
    }
    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)

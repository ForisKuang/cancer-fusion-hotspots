from cfh.algorithms import registry
from cfh.algorithms.confidence_stats import (
    ConfidenceStatsAlgorithm,
    default_confidence_stats_params,
    resolve_confidence_stats_params,
)
from cfh.genes.registry import load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


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


def test_real_pipeline_defaults_collapse_domain_states_and_use_tumor_variant_count():
    config = load_gene_config("RET")
    params = default_confidence_stats_params(config)

    # Regression guard: every event reaching the real pipeline has ALREADY
    # been filtered to Is_protein_fusion is True upstream (see
    # real_benchmark._is_target_protein_fusion), so an MLE default of
    # "Is_protein_fusion" is tautological -- guaranteed 100% success in
    # every non-empty group, always. It must never be the default outcome.
    assert params["outcome_field"] != "Is_protein_fusion"
    assert params["outcome_field"] == "Frame_status"
    assert params["outcome_success_value"] == "in-frame"

    statuses = ("retained", "retained", "lost", "lost", "disrupted", "disrupted")
    # Frame_status varies within each domain-retention group (including an
    # "unknown" that outcome_valid_values must exclude, not count as a
    # failure) precisely so the MLE below is a real, non-tautological
    # quantity, unlike the old Is_protein_fusion-based default.
    frame_statuses = ("in-frame", "out-of-frame", "in-frame", "out-of-frame", "in-frame", "unknown")
    events = [
        FusionEvent(
            Event_id=f"e{i}",
            Cohort="cohort1",
            Is_protein_fusion=True,
            Frame_status=frame_status,
            Tumor_variant_count=count,
        )
        for i, (count, frame_status) in enumerate(zip((30, 32, 4, 6, 8, 10), frame_statuses))
    ]
    features = [
        FusionFeature(
            Event_id=event.Event_id,
            Gene="RET",
            Domain_retention_flags={"kinase": status},
        )
        for event, status in zip(events, statuses, strict=True)
    ]

    result = ConfidenceStatsAlgorithm().run(events, features, config, params)

    assert result.Summary["group_a_label"] == "retained"
    assert result.Summary["group_b_label"] == "not_retained"
    assert result.Summary["n_group_a"] == 2
    assert result.Summary["n_group_b"] == 4
    assert result.Summary["mle"]["outcome_field"] == "Frame_status"
    # 1 of 2 retained events is in-frame -- a real proportion, not the
    # guaranteed-1.0 the old Is_protein_fusion default always produced.
    assert result.Summary["mle"]["groups"]["retained"]["n"] == 2
    assert result.Summary["mle"]["groups"]["retained"]["successes"] == 1
    assert result.Summary["mle"]["groups"]["retained"]["point_estimate"] != 1.0
    # The "unknown" not_retained event is excluded from the denominator
    # (3, not 4): it's undetermined, not a frame-status failure.
    assert result.Summary["mle"]["groups"]["not_retained"]["n"] == 3
    assert result.Summary["mle"]["groups"]["not_retained"]["successes"] == 2
    assert result.Summary["ttest"]["numeric_field"] == "Tumor_variant_count"
    assert result.Summary["ttest"]["p_value"] < 0.05


def test_default_params_can_be_overridden_by_real_pipeline_caller():
    custom = resolve_confidence_stats_params(
        load_gene_config("BRAF"),
        {"group_field": "Frame_status", "numeric_field": "Tumor_variant_count"},
    )

    assert custom["group_field"] == "Frame_status"
    assert "group_key" not in custom
    assert "group_value_map" not in custom
    assert "group_values" not in custom
    assert custom["numeric_field"] == "Tumor_variant_count"


def test_overriding_outcome_field_drops_incompatible_default_outcome_settings():
    """Overriding ``outcome_field`` away from the default ``Frame_status``
    must not silently keep the default's ``outcome_success_value``/
    ``outcome_valid_values`` -- those are meaningless (or actively
    misleading) for an unrelated field."""
    custom = resolve_confidence_stats_params(
        load_gene_config("BRAF"),
        {"outcome_field": "Is_antisense"},
    )

    assert custom["outcome_field"] == "Is_antisense"
    assert "outcome_success_value" not in custom
    assert "outcome_valid_values" not in custom


def test_outcome_success_value_required_for_a_non_boolean_field():
    """Without ``outcome_success_value``, a non-boolean ``outcome_field``
    (like a status string) falls back to ``bool(value)``, which is truthy
    for every non-empty value -- both "in-frame" and "out-of-frame" would
    register as "success". This is the bug the real default now avoids by
    always passing ``outcome_success_value`` explicitly."""
    events = [
        FusionEvent(Event_id="e0", Cohort="c", Frame_status="in-frame"),
        FusionEvent(Event_id="e1", Cohort="c", Frame_status="out-of-frame"),
        FusionEvent(Event_id="e2", Cohort="c", Frame_status="in-frame"),
        FusionEvent(Event_id="e3", Cohort="c", Frame_status="out-of-frame"),
    ]
    params_without_success_value = {
        "group_field": "Event_id",
        "group_values": ["e0", "e1"],
        "outcome_field": "Frame_status",
    }
    result = ConfidenceStatsAlgorithm().run(events, [], None, params_without_success_value)
    # Both "in-frame" and "out-of-frame" are truthy strings.
    assert result.Summary["mle"]["groups"]["e0"]["successes"] == 1
    assert result.Summary["mle"]["groups"]["e1"]["successes"] == 1

    params_with_success_value = {
        **params_without_success_value,
        "outcome_success_value": "in-frame",
    }
    result = ConfidenceStatsAlgorithm().run(events, [], None, params_with_success_value)
    assert result.Summary["mle"]["groups"]["e0"]["successes"] == 1
    assert result.Summary["mle"]["groups"]["e1"]["successes"] == 0


def test_outcome_valid_values_excludes_undetermined_observations_from_denominator():
    events = [
        FusionEvent(Event_id="e0", Cohort="c", Frame_status="in-frame"),
        FusionEvent(Event_id="e1", Cohort="c", Frame_status="unknown"),
        FusionEvent(Event_id="e2", Cohort="c", Frame_status="out-of-frame"),
    ]
    params = {
        "group_field": "Cohort",
        "group_values": ["c", "no-such-group"],
        "outcome_field": "Frame_status",
        "outcome_success_value": "in-frame",
        "outcome_valid_values": ["in-frame", "out-of-frame"],
    }

    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    group = result.Summary["mle"]["groups"]["c"]
    assert group["n"] == 2  # "unknown" excluded, not counted as a failure
    assert group["successes"] == 1


def test_mle_empty_group_emits_a_warning_matching_the_welch_skip_style():
    events = [
        FusionEvent(Event_id="e0", Cohort="c", Is_protein_fusion=True, Frame_status=None),
        FusionEvent(Event_id="e1", Cohort="other", Is_protein_fusion=True, Frame_status=None),
    ]
    params = {
        "group_field": "Cohort",
        "group_values": ["c", "other"],
        "outcome_field": "Frame_status",
        "outcome_success_value": "in-frame",
    }

    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert result.Summary["mle"]["groups"]["c"]["n"] == 0
    assert result.Summary["mle"]["groups"]["c"]["point_estimate"] is None
    assert any("MLE" in warning and "Frame_status" in warning for warning in result.Warnings), (
        result.Warnings
    )


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


def test_missing_group_field_is_a_clean_noop_not_a_raise():
    events = _make_events()

    result = ConfidenceStatsAlgorithm().run(events, [], None, {"outcome_field": "Is_antisense"})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "confidence_stats"
    assert result.Summary == {}
    assert result.Warnings
    assert "group_field" in result.Warnings[0]


def test_missing_outcome_and_numeric_field_is_a_clean_noop_not_a_raise():
    events = _make_events()

    result = ConfidenceStatsAlgorithm().run(events, [], None, {"group_field": "Is_protein_fusion"})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "confidence_stats"
    assert result.Summary == {}
    assert result.Warnings


def test_gene_config_with_no_group_field_configured_is_a_clean_noop():
    """A gene config lacking the required ``group_field`` param (e.g. a
    real full-cohort run for a gene that never opted into this optional
    comparison) must produce a no-op result, not an ``Algorithm failed``
    warning from a raised exception.
    """
    events = _make_events()
    config = load_gene_config("RET")

    result = ConfidenceStatsAlgorithm().run(events, [], config, {})

    assert result.Summary == {}
    assert result.Warnings == [
        "RET has no group_field configured; confidence-stats analysis was skipped."
    ]


def test_result_schema_matches_canonical_algorithm_result_fields():
    events = _make_events()
    params = {
        "group_field": "Is_protein_fusion",
        "outcome_field": "Is_antisense",
        "numeric_field": "Paired_end_read_support",
    }
    result = ConfidenceStatsAlgorithm().run(events, [], None, params)

    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)

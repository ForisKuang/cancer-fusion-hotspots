"""Generic statistical-confidence algorithm: binomial MLE/CI and Welch's t-test.

This algorithm is entirely gene-agnostic. It never references a specific
gene, domain, or biological concept -- it only knows how to (a) split
events/features into two groups based on a caller-specified field, (b)
optionally compute a binomial MLE point estimate and confidence interval for
a caller-specified boolean outcome within each group, and (c) optionally
compare a caller-specified numeric field between the two groups with
Welch's t-test. Domain-specific callers (e.g. the domain-retention module)
supply the field names via ``params``; this module has no knowledge of what
those fields mean biologically.

Opt-in via ``params["group_field"]`` (and at least one of
``params["outcome_field"]``/``params["numeric_field"]``): a gene run without
these configured produces a no-op result with no statistics computed, the
same graceful-skip pattern already used by ``exon_retention``/
``domain_disruption`` for genes that don't configure their respective
optional fields.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.mle import binomial_mle_confidence_interval
from cfh.stats.ttest import welch_t_test

ALGORITHM_NAME = "confidence_stats"
ALGORITHM_VERSION = "1.0.0"

_NO_SUCCESS_VALUE = object()
"""Sentinel distinguishing "no ``outcome_success_value`` given" from a
caller legitimately passing ``None`` as the success value."""


def default_confidence_stats_params(gene_config: GeneConfig) -> dict[str, Any]:
    """Return working, gene-agnostic defaults for a mapped fusion run.

    ``Domain_retention_flags`` is a mapping because one feature can carry
    several domain classifications.  The primary configured key merely
    selects which existing value to compare; no gene symbol or gene-specific
    statistical behavior is encoded here.  The algorithm is binary, so the
    biologically equivalent ``lost`` and ``disrupted`` states are explicitly
    collapsed into ``not_retained``.

    The MLE outcome defaults to ``Frame_status == "in-frame"`` (via
    ``outcome_success_value``), *not* ``Is_protein_fusion``: every event
    reaching the real benchmark has already been filtered to
    ``Is_protein_fusion is True`` upstream (see
    ``real_benchmark._is_target_protein_fusion``), so an MLE over
    ``Is_protein_fusion`` would be tautological -- guaranteed 100% success in
    every non-empty group, always, regardless of gene or cohort. That
    would feed a fixed, content-free CI width straight into
    ``composite_score``'s ``confidence_certainty`` sub-score, making larger
    groups look spuriously more "certain" for no biological reason.
    ``Frame_status`` is *not* filtered upstream and genuinely varies
    ("in-frame"/"out-of-frame"/"unknown"), so its per-group in-frame rate is
    a real, non-tautological quantity -- how often the retained-domain group
    is in-frame vs. how often the not-retained group is. ``unknown`` is
    excluded from the denominator via ``outcome_valid_values`` rather than
    counted as a failure, since it means "undetermined," not "out-of-frame".

    The numeric field defaults to ``Tumor_variant_count`` rather than
    ``Total_read_support``: on the live cBioPortal structural-variant API,
    ``tumorSplitReadCount``/``tumorPairedEndReadCount`` come back as the
    unavailable-value sentinel (``-1``, normalized to ``None``) for
    essentially every record, while ``tumorVariantCount`` is populated for
    the vast majority -- so it's the numeric field that actually lets the
    Welch comparison run on real data.
    """
    if not gene_config.key_domains:
        return {}
    primary_domain = gene_config.key_domains[0]
    domain_key = primary_domain.key or primary_domain.name
    return {
        "group_field": "Domain_retention_flags",
        "group_key": domain_key,
        "group_value_map": {
            "retained": "retained",
            "lost": "not_retained",
            "disrupted": "not_retained",
            "unknown": None,
        },
        "group_values": ["retained", "not_retained"],
        "outcome_field": "Frame_status",
        "outcome_success_value": "in-frame",
        "outcome_valid_values": ["in-frame", "out-of-frame"],
        "numeric_field": "Tumor_variant_count",
    }


def resolve_confidence_stats_params(
    gene_config: GeneConfig, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Merge caller overrides without retaining incompatible group/outcome defaults."""
    defaults = default_confidence_stats_params(gene_config)
    overrides = overrides or {}
    if overrides.get("group_field", defaults.get("group_field")) != defaults.get(
        "group_field"
    ):
        for key in ("group_key", "group_value_map", "group_values"):
            if key not in overrides:
                defaults.pop(key, None)
    if overrides.get("outcome_field", defaults.get("outcome_field")) != defaults.get(
        "outcome_field"
    ):
        for key in ("outcome_success_value", "outcome_valid_values"):
            if key not in overrides:
                defaults.pop(key, None)
    return {**defaults, **overrides}


def _build_rows(events: list[FusionEvent], features: list[FusionFeature]) -> list[dict]:
    """Join events and features on ``Event_id`` into flat field-name -> value dicts.

    One row is emitted per (event, feature) pair sharing an ``Event_id``. An
    event with no matching feature (or vice versa) still produces a row with
    whatever fields it has; missing fields simply resolve to ``None`` when
    looked up. This lets a caller name a field from either model without
    the algorithm needing to know which model it came from.
    """
    features_by_event: dict[str, list[dict]] = {}
    for feature in features:
        features_by_event.setdefault(feature.Event_id, []).append(feature.model_dump())

    rows: list[dict] = []
    seen_event_ids: set[str] = set()
    for event in events:
        event_dict = event.model_dump()
        matched = features_by_event.get(event.Event_id)
        if matched:
            for feature_dict in matched:
                rows.append({**event_dict, **feature_dict})
        else:
            rows.append(dict(event_dict))
        seen_event_ids.add(event.Event_id)

    for feature in features:
        if feature.Event_id not in seen_event_ids:
            rows.append(feature.model_dump())

    return rows


def _split_groups(
    rows: list[dict],
    group_field: str,
    group_values: Optional[list] = None,
    group_key: str | None = None,
    group_value_map: dict[Any, Any] | None = None,
) -> tuple[tuple[Any, list[dict]], tuple[Any, list[dict]]]:
    """Split rows into two groups by the distinct values of ``group_field``."""
    def group_value(row: dict) -> Any:
        value = row.get(group_field)
        if group_key is not None:
            value = value.get(group_key) if isinstance(value, dict) else None
        if group_value_map is not None:
            value = group_value_map.get(value, value)
        return value

    if group_values is not None:
        if len(group_values) != 2:
            raise ValueError("group_values must contain exactly two values")
        val_a, val_b = group_values
    else:
        distinct = sorted(
            {group_value(row) for row in rows if group_value(row) is not None},
            key=str,
        )
        if len(distinct) != 2:
            raise ValueError(
                f"group_field {group_field!r} must have exactly two distinct non-null "
                f"values across the input rows; found {distinct}"
            )
        val_a, val_b = distinct

    group_a = [row for row in rows if group_value(row) == val_a]
    group_b = [row for row in rows if group_value(row) == val_b]
    return (val_a, group_a), (val_b, group_b)


def _mle_block(
    group_a: tuple[Any, list[dict]],
    group_b: tuple[Any, list[dict]],
    outcome_field: str,
    confidence: float,
    method: str,
    outcome_success_value: Any = _NO_SUCCESS_VALUE,
    outcome_valid_values: Optional[list] = None,
) -> tuple[dict, list[str]]:
    """Compute the per-group binomial MLE/CI for ``outcome_field``.

    Success is ``value == outcome_success_value`` when that's given;
    otherwise falls back to ``bool(value)`` for a plain boolean field. A
    caller-facing string field (e.g. ``Frame_status``) is truthy for *every*
    non-empty value, so without ``outcome_success_value`` an MLE over such a
    field would silently compare the wrong thing rather than error --
    callers using a non-boolean ``outcome_field`` must pass it explicitly.

    ``outcome_valid_values``, when given, restricts the denominator to those
    values (e.g. excluding an "unknown"/undetermined state from counting as
    a failure) in addition to excluding ``None``.
    """
    label_a, rows_a = group_a
    label_b, rows_b = group_b
    groups = {}
    warnings: list[str] = []
    for label, rows in ((label_a, rows_a), (label_b, rows_b)):
        outcomes = [row.get(outcome_field) for row in rows if row.get(outcome_field) is not None]
        if outcome_valid_values is not None:
            outcomes = [value for value in outcomes if value in outcome_valid_values]
        n = len(outcomes)
        if outcome_success_value is not _NO_SUCCESS_VALUE:
            successes = sum(1 for value in outcomes if value == outcome_success_value)
        else:
            successes = sum(1 for value in outcomes if bool(value))
        if n:
            ci = binomial_mle_confidence_interval(
                successes, n, confidence=confidence, method=method
            )
        else:
            ci = {"point_estimate": None, "ci_low": None, "ci_high": None, "method": method}
            warnings.append(
                f"MLE binomial CI for {outcome_field!r} was not applicable for group "
                f"{label!r}; it had no valid non-null observations."
            )
        groups[str(label)] = {"n": n, "successes": successes, **ci}
    return {"outcome_field": outcome_field, "groups": groups}, warnings


def _ttest_block(
    group_a: tuple[Any, list[dict]],
    group_b: tuple[Any, list[dict]],
    numeric_field: str,
) -> dict:
    label_a, rows_a = group_a
    label_b, rows_b = group_b
    values_a = [row.get(numeric_field) for row in rows_a if row.get(numeric_field) is not None]
    values_b = [row.get(numeric_field) for row in rows_b if row.get(numeric_field) is not None]
    result = (
        welch_t_test(values_a, values_b)
        if len(values_a) >= 2 and len(values_b) >= 2
        else {
            "t_statistic": None,
            "p_value": None,
            "df": None,
            "mean_a": sum(values_a) / len(values_a) if values_a else None,
            "mean_b": sum(values_b) / len(values_b) if values_b else None,
        }
    )
    return {
        "numeric_field": numeric_field,
        "group_a_label": str(label_a),
        "group_b_label": str(label_b),
        "n_a": len(values_a),
        "n_b": len(values_b),
        **result,
    }


@register(ALGORITHM_NAME)
class ConfidenceStatsAlgorithm(Algorithm):
    """Gene-agnostic corroborating-confidence algorithm.

    Expected ``params`` keys (all optional -- a gene run without
    ``group_field`` set, or without at least one of ``outcome_field``/
    ``numeric_field``, gets a clean no-op result instead of this analysis):
        group_field (str): boolean/categorical field (on ``FusionEvent`` or
            ``FusionFeature``) defining the two groups being compared.
        group_values (list, optional): explicit [value_a, value_b] to use
            for the two groups; defaults to the two distinct non-null
            values observed for ``group_field``.
        group_key (str, optional): when ``group_field`` contains a mapping,
            compare the value at this key.
        group_value_map (dict, optional): map observed group values before
            splitting, for example to collapse several states into one.
        outcome_field (str, optional): field whose per-group success
            proportion is estimated via binomial MLE/CI. Omit to skip the
            MLE block entirely.
        outcome_success_value (any, optional): value that counts as
            "success" for ``outcome_field``. Required for a non-boolean
            ``outcome_field`` (e.g. a string like ``Frame_status``) --
            without it, every non-empty string is truthy and the MLE would
            silently compare the wrong thing. Omit only for a genuinely
            boolean ``outcome_field``, where it falls back to ``bool(value)``.
        outcome_valid_values (list, optional): restrict the MLE denominator
            to these ``outcome_field`` values, excluding e.g. an "unknown"/
            undetermined state from counting as a failure.
        confidence (float, optional): confidence level for the MLE
            interval, default 0.95.
        mle_method (str, optional): ``"wilson"`` (default) or
            ``"clopper_pearson"``.
        numeric_field (str, optional): numeric field compared between the
            two groups with Welch's t-test. Omit to skip the t-test block
            entirely.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: Optional[GeneConfig],
        params: dict,
    ) -> AlgorithmResult:
        group_field = params.get("group_field")
        outcome_field = params.get("outcome_field")
        numeric_field = params.get("numeric_field")
        gene_name = (
            (gene_config.gene_symbol or gene_config.gene_pair) if gene_config else None
        ) or "This gene"
        if not group_field:
            return self._no_op_result(
                f"{gene_name} has no group_field configured; confidence-stats "
                "analysis was skipped."
            )
        if not outcome_field and not numeric_field:
            return self._no_op_result(
                f"{gene_name} configured a group_field but no outcome_field or "
                "numeric_field; confidence-stats analysis was skipped."
            )

        confidence = params.get("confidence", 0.95)
        mle_method = params.get("mle_method", "wilson")
        group_values = params.get("group_values")
        group_key = params.get("group_key")
        group_value_map = params.get("group_value_map")
        outcome_success_value = params.get("outcome_success_value", _NO_SUCCESS_VALUE)
        outcome_valid_values = params.get("outcome_valid_values")

        rows = _build_rows(events, features)
        group_a, group_b = _split_groups(
            rows, group_field, group_values, group_key, group_value_map
        )

        summary: dict = {
            "group_field": group_field,
            "group_a_label": str(group_a[0]),
            "group_b_label": str(group_b[0]),
            "n_group_a": len(group_a[1]),
            "n_group_b": len(group_b[1]),
        }

        warnings: list[str] = []
        if outcome_field:
            summary["mle"], mle_warnings = _mle_block(
                group_a,
                group_b,
                outcome_field,
                confidence,
                mle_method,
                outcome_success_value,
                outcome_valid_values,
            )
            warnings.extend(mle_warnings)

        if numeric_field:
            summary["ttest"] = _ttest_block(group_a, group_b, numeric_field)
            if summary["ttest"]["p_value"] is None:
                warnings.append(
                    f"Welch's t-test for {numeric_field!r} was not applicable; "
                    "each group requires at least two non-null observations."
                )

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters=params,
            Summary=summary,
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )

    @staticmethod
    def _no_op_result(warning: str) -> AlgorithmResult:
        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={},
            Summary={},
            Warnings=[warning],
            Created_at=datetime.now(timezone.utc),
        )

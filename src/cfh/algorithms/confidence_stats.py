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
) -> tuple[tuple[Any, list[dict]], tuple[Any, list[dict]]]:
    """Split rows into two groups by the distinct values of ``group_field``."""
    if group_values is not None:
        if len(group_values) != 2:
            raise ValueError("group_values must contain exactly two values")
        val_a, val_b = group_values
    else:
        distinct = sorted(
            {row.get(group_field) for row in rows if row.get(group_field) is not None},
            key=str,
        )
        if len(distinct) != 2:
            raise ValueError(
                f"group_field {group_field!r} must have exactly two distinct non-null "
                f"values across the input rows; found {distinct}"
            )
        val_a, val_b = distinct

    group_a = [row for row in rows if row.get(group_field) == val_a]
    group_b = [row for row in rows if row.get(group_field) == val_b]
    return (val_a, group_a), (val_b, group_b)


def _mle_block(
    group_a: tuple[Any, list[dict]],
    group_b: tuple[Any, list[dict]],
    outcome_field: str,
    confidence: float,
    method: str,
) -> dict:
    label_a, rows_a = group_a
    label_b, rows_b = group_b
    groups = {}
    for label, rows in ((label_a, rows_a), (label_b, rows_b)):
        outcomes = [row.get(outcome_field) for row in rows if row.get(outcome_field) is not None]
        n = len(outcomes)
        successes = sum(1 for value in outcomes if bool(value))
        ci = binomial_mle_confidence_interval(successes, n, confidence=confidence, method=method)
        groups[str(label)] = {"n": n, "successes": successes, **ci}
    return {"outcome_field": outcome_field, "groups": groups}


def _ttest_block(
    group_a: tuple[Any, list[dict]],
    group_b: tuple[Any, list[dict]],
    numeric_field: str,
) -> dict:
    label_a, rows_a = group_a
    label_b, rows_b = group_b
    values_a = [row.get(numeric_field) for row in rows_a if row.get(numeric_field) is not None]
    values_b = [row.get(numeric_field) for row in rows_b if row.get(numeric_field) is not None]
    result = welch_t_test(values_a, values_b)
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
        outcome_field (str, optional): boolean field whose per-group
            proportion is estimated via binomial MLE/CI. Omit to skip the
            MLE block entirely.
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
                f"{gene_name} has no group_field configured; confidence-stats analysis was skipped."
            )
        if not outcome_field and not numeric_field:
            return self._no_op_result(
                f"{gene_name} configured a group_field but no outcome_field or "
                "numeric_field; confidence-stats analysis was skipped."
            )

        confidence = params.get("confidence", 0.95)
        mle_method = params.get("mle_method", "wilson")
        group_values = params.get("group_values")

        rows = _build_rows(events, features)
        group_a, group_b = _split_groups(rows, group_field, group_values)

        summary: dict = {
            "group_field": group_field,
            "group_a_label": str(group_a[0]),
            "group_b_label": str(group_b[0]),
            "n_group_a": len(group_a[1]),
            "n_group_b": len(group_b[1]),
        }

        if outcome_field:
            summary["mle"] = _mle_block(group_a, group_b, outcome_field, confidence, mle_method)

        if numeric_field:
            summary["ttest"] = _ttest_block(group_a, group_b, numeric_field)

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters=params,
            Summary=summary,
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

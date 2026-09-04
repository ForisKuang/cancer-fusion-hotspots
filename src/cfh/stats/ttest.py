"""Generic Welch's t-test utility for comparing two independent numeric samples.

Gene-agnostic: callers pass two lists of numeric observations for whatever
groups they are comparing (e.g. read-support values split by a retention
flag).
"""

from __future__ import annotations

import math
import statistics

from scipy.stats import ttest_ind


def _welch_satterthwaite_df(group_a: list[float], group_b: list[float]) -> float:
    """Compute the Welch-Satterthwaite degrees of freedom independently of scipy."""
    n_a, n_b = len(group_a), len(group_b)
    var_a = statistics.variance(group_a)
    var_b = statistics.variance(group_b)
    se_a2 = var_a / n_a
    se_b2 = var_b / n_b
    numerator = (se_a2 + se_b2) ** 2
    denominator = (se_a2**2) / (n_a - 1) + (se_b2**2) / (n_b - 1)
    if denominator == 0:
        # Both groups have zero variance -- the Satterthwaite formula is 0/0.
        # scipy.stats.ttest_ind falls back to df=1.0 in this case; match it.
        return 1.0
    return numerator / denominator


def welch_t_test(group_a: list[float], group_b: list[float]) -> dict:
    """Run Welch's t-test (unequal variances) comparing two independent samples.

    Uses ``scipy.stats.ttest_ind(..., equal_var=False)`` for the t-statistic
    and p-value, but computes the Welch-Satterthwaite degrees of freedom
    explicitly rather than trusting scipy's internals, and cross-checks the
    two against each other when scipy exposes its own ``df`` attribute.

    Args:
        group_a: numeric observations for the first group (n >= 2).
        group_b: numeric observations for the second group (n >= 2).

    Returns:
        dict with keys ``t_statistic``, ``p_value``, ``df``, ``mean_a``, ``mean_b``.
    """
    if len(group_a) < 2 or len(group_b) < 2:
        raise ValueError("each group must have at least 2 observations")

    result = ttest_ind(group_a, group_b, equal_var=False)
    df = _welch_satterthwaite_df(group_a, group_b)

    scipy_df = getattr(result, "df", None)
    if scipy_df is not None and not math.isclose(df, float(scipy_df), rel_tol=1e-9, abs_tol=1e-9):
        raise AssertionError(
            f"Independently computed Welch-Satterthwaite df ({df}) does not match "
            f"scipy's reported df ({scipy_df})"
        )

    return {
        "t_statistic": float(result.statistic),
        "p_value": float(result.pvalue),
        "df": df,
        "mean_a": statistics.mean(group_a),
        "mean_b": statistics.mean(group_b),
    }

import math

import pytest

from cfh.stats.mle import binomial_mle_confidence_interval


def _independent_wilson_interval(successes: int, n: int, confidence: float) -> tuple[float, float]:
    """Wilson score interval, computed from scratch (not copy-pasted from the
    implementation under test), to independently verify the module's output."""
    from scipy.stats import norm

    p_hat = successes / n
    z = norm.ppf(1 - (1 - confidence) / 2)
    denom = 1 + (z**2) / n
    center = (p_hat + (z**2) / (2 * n)) / denom
    margin = (z / denom) * math.sqrt(p_hat * (1 - p_hat) / n + (z**2) / (4 * n**2))
    return max(0.0, center - margin), min(1.0, center + margin)


def test_wilson_interval_matches_independent_formula_and_brackets_point_estimate():
    result = binomial_mle_confidence_interval(30, 33, method="wilson")

    assert result["point_estimate"] == 30 / 33
    assert result["method"] == "wilson"

    expected_low, expected_high = _independent_wilson_interval(30, 33, 0.95)
    assert result["ci_low"] == pytest.approx(expected_low)
    assert result["ci_high"] == pytest.approx(expected_high)

    assert 0.0 <= result["ci_low"] <= result["point_estimate"] <= result["ci_high"] <= 1.0


@pytest.mark.parametrize("successes,n", [(0, 10), (10, 10)])
def test_clopper_pearson_edge_cases_return_finite_valid_bounds(successes, n):
    result = binomial_mle_confidence_interval(successes, n, method="clopper_pearson")

    assert result["method"] == "clopper_pearson"
    assert result["point_estimate"] == successes / n
    assert math.isfinite(result["ci_low"])
    assert math.isfinite(result["ci_high"])
    assert 0.0 <= result["ci_low"] <= result["ci_high"] <= 1.0


def test_mle_point_estimate_is_successes_over_n():
    result = binomial_mle_confidence_interval(7, 20, method="wilson")
    assert result["point_estimate"] == 7 / 20


def test_unknown_method_raises():
    with pytest.raises(ValueError):
        binomial_mle_confidence_interval(1, 2, method="bogus")


def test_invalid_counts_raise():
    with pytest.raises(ValueError):
        binomial_mle_confidence_interval(-1, 10)
    with pytest.raises(ValueError):
        binomial_mle_confidence_interval(11, 10)
    with pytest.raises(ValueError):
        binomial_mle_confidence_interval(1, 0)

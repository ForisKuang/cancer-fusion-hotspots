"""Generic maximum-likelihood confidence-interval utilities for binomial proportions.

Gene-agnostic: callers pass raw success/trial counts for whatever boolean
outcome they are studying (e.g. "was a given domain retained?").
"""

from __future__ import annotations

import math

from scipy.stats import beta as beta_dist
from scipy.stats import norm


def _wilson_score_interval(successes: int, n: int, confidence: float) -> tuple[float, float]:
    """Closed-form Wilson score interval for a binomial proportion."""
    p_hat = successes / n
    z = norm.ppf(1 - (1 - confidence) / 2)
    z2 = z * z
    denom = 1 + z2 / n
    center = (p_hat + z2 / (2 * n)) / denom
    half_width = (z / denom) * math.sqrt((p_hat * (1 - p_hat) / n) + (z2 / (4 * n * n)))
    ci_low = max(0.0, center - half_width)
    ci_high = min(1.0, center + half_width)
    return ci_low, ci_high


def _clopper_pearson_interval(successes: int, n: int, confidence: float) -> tuple[float, float]:
    """Exact Clopper-Pearson interval for a binomial proportion, via the beta distribution."""
    alpha = 1 - confidence
    if successes == 0:
        ci_low = 0.0
    else:
        ci_low = beta_dist.ppf(alpha / 2, successes, n - successes + 1)
    if successes == n:
        ci_high = 1.0
    else:
        ci_high = beta_dist.ppf(1 - alpha / 2, successes + 1, n - successes)
    return float(ci_low), float(ci_high)


def binomial_mle_confidence_interval(
    successes: int,
    n: int,
    confidence: float = 0.95,
    method: str = "wilson",
) -> dict:
    """Compute a maximum-likelihood point estimate and confidence interval for a
    binomial proportion.

    The maximum-likelihood estimator (MLE) of a binomial proportion is simply
    ``successes / n`` -- this is the value that maximizes the binomial
    likelihood function over the observed data, and is computed explicitly
    below rather than being treated as an incidental "rate".

    Args:
        successes: number of successes observed (0 <= successes <= n).
        n: number of trials (n > 0).
        confidence: confidence level in (0, 1), default 0.95.
        method: ``"wilson"`` (closed-form Wilson score interval, default) or
            ``"clopper_pearson"`` (exact interval via the beta distribution).

    Returns:
        dict with keys ``point_estimate``, ``ci_low``, ``ci_high``, ``method``.
    """
    if n <= 0:
        raise ValueError("n must be a positive integer")
    if not (0 <= successes <= n):
        raise ValueError("successes must satisfy 0 <= successes <= n")
    if not (0 < confidence < 1):
        raise ValueError("confidence must be in (0, 1)")

    # Maximum-likelihood estimate of the binomial success probability.
    point_estimate = successes / n

    if method == "wilson":
        ci_low, ci_high = _wilson_score_interval(successes, n, confidence)
    elif method == "clopper_pearson":
        ci_low, ci_high = _clopper_pearson_interval(successes, n, confidence)
    else:
        raise ValueError(f"Unknown method {method!r}; expected 'wilson' or 'clopper_pearson'")

    # Floating-point residue can otherwise place the interval bound on the
    # wrong side of the point estimate at the successes=0 or successes=n
    # boundary (e.g. ci_high == 0.9999999999999999 when point_estimate == 1.0).
    ci_low = min(ci_low, point_estimate)
    ci_high = max(ci_high, point_estimate)

    return {
        "point_estimate": point_estimate,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "method": method,
    }

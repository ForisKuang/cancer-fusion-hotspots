"""Adaptive permutation-count budgeting shared by every algorithm that runs
a label/breakpoint permutation test (``domain_retention``,
``domain_disruption``, ``cutpoint_detection``).

Running the full configured ``n_permutations`` (thousands) for every gene
in a genome-wide scan is wasteful when a result is nowhere near the
significance threshold. Adaptive mode instead runs a small permutation
budget first and only escalates to the full budget when that small-N
result is "borderline" -- close enough to the significance threshold that
the small-N Monte Carlo estimate's own noise could plausibly flip the
significance call.

Both the small and full runs reuse the exact same ``seed``, so a given
``(seed, n_permutations)`` pair always reproduces the same result -- exactly
as deterministic as the existing non-adaptive path. Escalating just means
calling the same seeded permutation routine again with a larger sample
count, never continuing an unseeded stream.

Adaptive mode is strictly opt-in (``params["adaptive"]``, default
``False``): every existing caller that never sets these keys keeps its
exact current, already-tested behavior and output shape -- this module only
adds a new code path alongside it.
"""

from __future__ import annotations

DEFAULT_SMALL_N_PERMUTATIONS = 100
DEFAULT_SIGNIFICANCE_THRESHOLD = 0.05
DEFAULT_BORDERLINE_FACTOR = 2.0


def is_borderline(
    p_value: float | None,
    *,
    threshold: float = DEFAULT_SIGNIFICANCE_THRESHOLD,
    factor: float = DEFAULT_BORDERLINE_FACTOR,
) -> bool:
    """Whether ``p_value`` sits within ``factor``x of ``threshold`` on either side.

    ``None`` (no p-value could be computed, e.g. an indeterminable cutpoint
    scan) is never borderline -- there is nothing to escalate.
    """
    if p_value is None:
        return False
    if factor <= 0:
        raise ValueError(f"factor must be positive; got {factor!r}")
    if not (0.0 <= threshold <= 1.0):
        raise ValueError(f"threshold must be between 0 and 1; got {threshold!r}")
    lower = threshold / factor
    upper = min(1.0, threshold * factor)
    return lower <= p_value <= upper


def resolve_permutation_budget(params: dict, *, default_full_n: int) -> dict:
    """Resolve the adaptive-permutation knobs from an algorithm's ``params``.

    Returns a dict with keys ``adaptive``, ``small_n``, ``full_n``,
    ``threshold``, ``factor``. Every key is optional in ``params`` and
    defaults to non-adaptive (``adaptive=False``) with ``full_n`` equal to
    whatever ``params["n_permutations"]`` already resolves to (or
    ``default_full_n``), so an algorithm's existing non-adaptive callers see
    no change at all.

    Recognized ``params`` keys:
        adaptive (bool): opt in to adaptive budgeting. Default ``False``.
        n_permutations (int): the full/escalated permutation count -- the
            same key non-adaptive callers already use.
        n_permutations_small (int): the small initial permutation count.
            Default :data:`DEFAULT_SMALL_N_PERMUTATIONS`.
        significance_threshold (float): the significance level used to
            judge "borderline". Default :data:`DEFAULT_SIGNIFICANCE_THRESHOLD`.
        borderline_factor (float): how many multiples of the threshold
            around it still count as borderline. Default
            :data:`DEFAULT_BORDERLINE_FACTOR`.
    """
    return {
        "adaptive": bool(params.get("adaptive", False)),
        "small_n": int(params.get("n_permutations_small", DEFAULT_SMALL_N_PERMUTATIONS)),
        "full_n": int(params.get("n_permutations", default_full_n)),
        "threshold": float(params.get("significance_threshold", DEFAULT_SIGNIFICANCE_THRESHOLD)),
        "factor": float(params.get("borderline_factor", DEFAULT_BORDERLINE_FACTOR)),
    }

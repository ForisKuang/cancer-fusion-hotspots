"""Multiple-testing corrections for cross-gene analyses."""

from __future__ import annotations

import math
from collections.abc import Sequence

Hypothesis = tuple[str, str, float]
AdjustedHypothesis = tuple[str, str, float, float]


def benjamini_hochberg(hypotheses: Sequence[Hypothesis]) -> list[AdjustedHypothesis]:
    """Apply Benjamini-Hochberg FDR correction, preserving input order.

    Each input is ``(gene, algorithm, p_value)``. The returned tuples append
    the corresponding BH-adjusted q-value. All supplied hypotheses form one
    correction family.
    """
    validated: list[Hypothesis] = []
    for gene, algorithm, p_value in hypotheses:
        if isinstance(p_value, bool):
            raise ValueError(
                f"p-value for {gene}/{algorithm} must be finite and between 0 and 1; "
                f"got {p_value!r}"
            )
        value = float(p_value)
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(
                f"p-value for {gene}/{algorithm} must be finite and between 0 and 1; "
                f"got {p_value!r}"
            )
        validated.append((gene, algorithm, value))

    count = len(validated)
    if count == 0:
        return []

    sorted_indices = sorted(range(count), key=lambda index: validated[index][2])
    adjusted = [0.0] * count
    running_minimum = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = sorted_indices[rank_index]
        rank = rank_index + 1
        candidate = validated[original_index][2] * count / rank
        running_minimum = min(running_minimum, candidate)
        adjusted[original_index] = running_minimum

    return [
        (gene, algorithm, p_value, adjusted[index])
        for index, (gene, algorithm, p_value) in enumerate(validated)
    ]

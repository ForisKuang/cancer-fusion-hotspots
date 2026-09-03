import math

import pytest

from cfh.stats.cutpoint_scan import (
    candidate_cutpoints,
    detect_cutpoint,
    fishers_cutpoint_test,
    scan_cutpoints,
)


def _clean_separation():
    """Retained breakpoints cluster below 120, lost breakpoints above 195 --
    an unambiguous single boundary between them."""
    positions = [100, 105, 110, 115, 120, 195, 200, 205, 210, 215]
    statuses = ["retained"] * 5 + ["lost"] * 5
    return positions, statuses


def test_candidate_cutpoints_excludes_the_largest_observed_position():
    assert candidate_cutpoints([300, 100, 200, 100]) == [100, 200]


def test_scan_cutpoints_builds_correct_contingency_counts():
    positions, statuses = _clean_separation()
    scan = scan_cutpoints(positions, statuses)

    perfect_row = next(row for row in scan if row["cutpoint"] == 120)
    assert perfect_row["n_positive_at_or_below"] == 5
    assert perfect_row["n_positive_above"] == 0
    assert perfect_row["n_negative_at_or_below"] == 0
    assert perfect_row["n_negative_above"] == 5
    assert perfect_row["odds_ratio"] == float("inf")

    weaker_row = next(row for row in scan if row["cutpoint"] == 115)
    assert weaker_row["n_positive_at_or_below"] == 4
    assert weaker_row["n_positive_above"] == 1
    assert weaker_row["neg_log10_p_value"] < perfect_row["neg_log10_p_value"]


def test_fishers_cutpoint_test_matches_scipy_two_sided():
    from scipy.stats import fisher_exact

    table = [[6, 0], [0, 8]]
    odds_ratio, p_value = fishers_cutpoint_test(table)
    expected_odds_ratio, expected_p_value = fisher_exact(table, alternative="two-sided")
    assert odds_ratio == float(expected_odds_ratio)
    assert p_value == pytest.approx(float(expected_p_value))


def test_fishers_cutpoint_test_rejects_bad_tables():
    with pytest.raises(ValueError):
        fishers_cutpoint_test([[1, 2, 3], [4, 5, 6]])
    with pytest.raises(ValueError):
        fishers_cutpoint_test([[-1, 2], [3, 4]])


def test_clean_single_boundary_recovers_the_true_cutpoint():
    positions, statuses = _clean_separation()

    result = detect_cutpoint(positions, statuses, seed=42, n_permutations=500)

    assert result["determinable"] is True
    assert result["reason"] is None
    assert result["best_cutpoint"] == 120
    assert result["observed_odds_ratio"] == float("inf")
    assert math.isfinite(result["observed_statistic"]) or result["observed_statistic"] == math.inf
    assert 0.0 <= result["corrected_p_value"] <= 1.0
    assert len(result["scan"]) == len(candidate_cutpoints(positions))
    assert len(result["null_max_statistics"]) == 500


def test_permutation_correction_is_never_more_optimistic_than_the_naive_scan_minimum():
    """Sanity check on the multiple-comparisons correction: since the
    observed statistic is the *best of many* candidate cutpoints, its
    empirical (permutation-corrected) p-value should never be smaller than
    the naive per-cutpoint Fisher p-value at that same cutpoint -- and, with
    a small permutation budget, is bounded below by 1/(B+1)."""
    positions, statuses = _clean_separation()
    n_permutations = 200

    result = detect_cutpoint(positions, statuses, seed=7, n_permutations=n_permutations)

    assert result["corrected_p_value"] >= result["observed_p_value"]
    assert result["corrected_p_value"] >= 1 / (n_permutations + 1)
    # A real, strong signal should still come out as clearly significant
    # after correction.
    assert result["corrected_p_value"] < 0.1


def test_permutation_is_bit_identical_for_a_fixed_seed():
    positions, statuses = _clean_separation()
    first = detect_cutpoint(positions, statuses, seed=1729, n_permutations=250)
    second = detect_cutpoint(positions, statuses, seed=1729, n_permutations=250)
    assert first == second


def test_no_signal_case_does_not_crash_and_stays_bounded():
    """Breakpoints and outcome are interleaved with no positional pattern --
    scanning must still complete and return a valid, bounded result."""
    positions = list(range(100, 120))
    statuses = ["retained" if position % 2 == 0 else "lost" for position in positions]

    result = detect_cutpoint(positions, statuses, seed=3, n_permutations=200)

    assert result["determinable"] is True
    assert 0.0 <= result["corrected_p_value"] <= 1.0
    assert result["best_cutpoint"] in positions


@pytest.mark.parametrize(
    "positions,statuses,expected_reason_snippet",
    [
        ([100, 110, 120], ["retained", "lost", "retained"], "fewer than 4"),
        ([100, 100, 100, 100], ["retained", "retained", "lost", "lost"], "distinct breakpoint"),
        ([100, 110, 120, 130], ["retained"] * 4, "single outcome class"),
    ],
)
def test_degenerate_inputs_return_not_determinable_without_raising(
    positions, statuses, expected_reason_snippet
):
    result = detect_cutpoint(positions, statuses, seed=42, n_permutations=100)

    assert result["determinable"] is False
    assert expected_reason_snippet in result["reason"]
    assert result["best_cutpoint"] is None
    assert result["corrected_p_value"] is None
    assert result["scan"] == [] or all(isinstance(row, dict) for row in result["scan"])


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        detect_cutpoint([100, 200], ["retained"], seed=1, n_permutations=10)


def test_non_positive_permutations_raise():
    positions = [100, 110, 120, 130]
    statuses = ["retained", "retained", "lost", "lost"]
    with pytest.raises(ValueError):
        detect_cutpoint(positions, statuses, n_permutations=0)

import random

import pytest
from scipy.stats import ttest_ind

from cfh.stats.ttest import welch_t_test


def test_p_value_matches_scipy_ttest_ind_independently_recomputed():
    group_a = [1, 2, 3, 4, 5]
    group_b = [10, 20, 30, 40, 50]

    result = welch_t_test(group_a, group_b)

    expected = ttest_ind(group_a, group_b, equal_var=False)
    assert result["p_value"] == pytest.approx(expected.pvalue)
    assert result["t_statistic"] == pytest.approx(expected.statistic)
    assert result["mean_a"] == pytest.approx(3.0)
    assert result["mean_b"] == pytest.approx(30.0)


def test_reported_df_matches_scipy_welch_satterthwaite_df():
    group_a = [1, 2, 3, 4, 5, 6]
    group_b = [10, 20, 30, 40, 50]

    result = welch_t_test(group_a, group_b)
    expected = ttest_ind(group_a, group_b, equal_var=False)

    assert hasattr(expected, "df")
    assert result["df"] == pytest.approx(float(expected.df))


def test_positive_and_negative_controls():
    # Positive control: same variance, very different means -> small p-value.
    high_signal_a = [10.0, 10.5, 9.5, 10.2, 9.8, 10.1]
    high_signal_b = [50.0, 50.5, 49.5, 50.2, 49.8, 50.1]
    signal_result = welch_t_test(high_signal_a, high_signal_b)
    assert signal_result["p_value"] < 0.05

    # Negative control: same underlying distribution, just reordered/resampled.
    rng = random.Random(42)
    base = [10.0, 12.0, 9.0, 11.0, 10.5, 9.5, 11.5, 10.2, 9.8, 10.8]
    group_a = rng.sample(base, k=len(base))
    group_b = rng.sample(base, k=len(base))
    noise_result = welch_t_test(group_a, group_b)
    assert noise_result["p_value"] > 0.05


def test_requires_at_least_two_observations_per_group():
    with pytest.raises(ValueError):
        welch_t_test([1.0], [1.0, 2.0])

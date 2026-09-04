import pytest

from cfh.stats.multiple_testing import benjamini_hochberg


def test_benjamini_hochberg_matches_hand_computed_textbook_example():
    hypotheses = [
        ("GENE_A", "method", 0.01),
        ("GENE_B", "method", 0.04),
        ("GENE_C", "method", 0.03),
        ("GENE_D", "method", 0.002),
    ]

    adjusted = benjamini_hochberg(hypotheses)

    # Sorted p-values .002, .01, .03, .04 have monotone BH q-values
    # .008, .02, .04, .04. The function restores original input order.
    assert [row[3] for row in adjusted] == pytest.approx([0.02, 0.04, 0.04, 0.008])
    assert [row[:3] for row in adjusted] == hypotheses


def test_benjamini_hochberg_handles_empty_input_and_tied_p_values():
    assert benjamini_hochberg([]) == []
    adjusted = benjamini_hochberg(
        [("GENE_A", "method", 0.5), ("GENE_B", "method", 0.5)]
    )
    assert [row[3] for row in adjusted] == [0.5, 0.5]


@pytest.mark.parametrize("p_value", [-0.01, 1.01, float("nan"), float("inf"), True])
def test_benjamini_hochberg_rejects_invalid_p_values(p_value):
    with pytest.raises(ValueError, match="between 0 and 1"):
        benjamini_hochberg([("GENE_A", "method", p_value)])

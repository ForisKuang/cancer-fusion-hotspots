"""Unit tests for the gene-agnostic genome-wide Manhattan/volcano SVG
(:mod:`cfh.reporting.manhattan`), plus a validation against the real,
already-committed ``msk_impact_50k_2026`` cohort-scan run.

The real-data test is the one that matters most: it checks the SVG's own
generated coordinates (not just an eyeballed picture) to confirm the
already-verified finding from the locus-validated rerun of that cohort scan
-- ETV6 (q~0.0043) is FDR-significant and plotted above the dashed q=0.05
threshold line, while RET (q~0.1198) and BRAF (q~0.2356) are not significant
and plotted below it. (RET's significance call flipped after the locus
validation fix corrected its breakpoint mapping -- see the ``real_benchmark``
rerun that produced this run.)
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pytest

from cfh.reporting.manhattan import render_manhattan_svg

REPO_ROOT = Path(__file__).parent.parent
REAL_COHORT_SCAN_SUMMARY_JSON = (
    REPO_ROOT
    / "runs"
    / "cohort-scan_msk_impact_50k_2026_20260904T005352Z"
    / "cohort_scan"
    / "summary.json"
)


def _row(
    gene_symbol: str,
    *,
    q_value: float | None,
    composite_score: float | None,
    significant: bool = False,
) -> dict:
    return {
        "gene_symbol": gene_symbol,
        "min_fdr_adjusted_q_value": q_value,
        "top_composite_score": composite_score,
        "fdr_significant": significant,
    }


def _circle_cy(svg: str, gene_symbol: str) -> float:
    match = re.search(
        rf'<circle cx="([0-9.]+)" cy="([0-9.]+)"[^>]*data-gene="{re.escape(gene_symbol)}"', svg
    )
    assert match, f"no plotted circle found for {gene_symbol}"
    return float(match.group(2))


def _threshold_line_y(svg: str) -> float:
    match = re.search(r'id="fdr-threshold-line"[^>]*y1="([0-9.]+)"', svg)
    assert match, "no FDR threshold line found in SVG"
    return float(match.group(1))


def test_svg_is_well_formed_and_gene_agnostic():
    rows = [
        _row("GENEA", q_value=0.001, composite_score=0.9, significant=True),
        _row("GENEB", q_value=0.2, composite_score=0.5, significant=False),
        _row("GENEC", q_value=None, composite_score=None),  # never FDR-tested: excluded
    ]
    svg = render_manhattan_svg(rows, significance_level=0.05)

    assert svg.startswith('<svg xmlns="http://www.w3.org/2000/svg"')
    assert svg.endswith("</svg>")
    assert svg.count("<circle") - 2 == 2  # 2 plotted genes + 2 legend swatches
    assert 'data-gene="GENEA"' in svg
    assert 'data-gene="GENEB"' in svg
    assert 'data-gene="GENEC"' not in svg  # no q-value -> not plotted


def test_significant_points_are_above_the_dashed_threshold_line():
    rows = [
        _row("SIGGENE", q_value=0.001, composite_score=0.8, significant=True),
        _row("NOTSIGGENE", q_value=0.5, composite_score=0.2, significant=False),
    ]
    svg = render_manhattan_svg(rows, significance_level=0.05)

    threshold_y = _threshold_line_y(svg)
    assert _circle_cy(svg, "SIGGENE") < threshold_y  # smaller y = higher = above the line
    assert _circle_cy(svg, "NOTSIGGENE") > threshold_y


def test_significant_and_not_significant_points_use_different_colors():
    rows = [
        _row("SIGGENE", q_value=0.001, composite_score=0.8, significant=True),
        _row("NOTSIGGENE", q_value=0.5, composite_score=0.2, significant=False),
    ]
    svg = render_manhattan_svg(rows, significance_level=0.05)

    sig_circle = re.search(r'<circle[^>]*data-gene="SIGGENE"[^>]*/>', svg).group(0)
    not_sig_circle = re.search(r'<circle[^>]*data-gene="NOTSIGGENE"[^>]*/>', svg).group(0)
    sig_fill = re.search(r'fill="(#[0-9a-fA-F]+)"', sig_circle).group(1)
    not_sig_fill = re.search(r'fill="(#[0-9a-fA-F]+)"', not_sig_circle).group(1)
    assert sig_fill != not_sig_fill
    assert 'data-significant="true"' in sig_circle
    assert 'data-significant="false"' in not_sig_circle


def test_x_axis_ranks_by_descending_composite_score():
    rows = [
        _row("LOWSCORE", q_value=0.3, composite_score=0.1),
        _row("HIGHSCORE", q_value=0.3, composite_score=0.9),
        _row("MIDSCORE", q_value=0.3, composite_score=0.5),
    ]
    svg = render_manhattan_svg(rows, significance_level=0.05)

    def cx(gene_symbol: str) -> float:
        match = re.search(
            rf'<circle cx="([0-9.]+)"[^>]*data-gene="{gene_symbol}"',
            svg,
        )
        return float(match.group(1))

    assert cx("HIGHSCORE") < cx("MIDSCORE") < cx("LOWSCORE")


def test_labels_include_every_significant_gene_even_beyond_max_labels():
    rows = [
        _row(f"SIG{i}", q_value=0.001, composite_score=0.9 - i * 0.01, significant=True)
        for i in range(5)
    ] + [_row("NOTSIG", q_value=0.5, composite_score=0.1, significant=False)]
    svg = render_manhattan_svg(rows, significance_level=0.05, max_labels=3)

    for i in range(5):
        assert f">SIG{i}<" in svg


def test_labels_fill_remaining_slots_by_significance_not_composite_rank():
    rows = [
        _row("SIGGENE", q_value=0.001, composite_score=0.1, significant=True),
        _row("BESTQNOTSIG", q_value=0.06, composite_score=0.99, significant=False),
        _row("WORSTQNOTSIG", q_value=0.9, composite_score=0.5, significant=False),
    ]
    svg = render_manhattan_svg(rows, significance_level=0.05, max_labels=2)

    assert ">SIGGENE<" in svg
    assert ">BESTQNOTSIG<" in svg
    assert ">WORSTQNOTSIG<" not in svg


def test_empty_and_all_untested_rows_render_a_valid_placeholder_svg():
    assert render_manhattan_svg([]).endswith("</svg>")
    only_untested = [_row("GENEA", q_value=None, composite_score=None)]
    svg = render_manhattan_svg(only_untested)
    assert "No scanned gene produced an FDR-adjusted q-value" in svg
    assert "data-gene=" not in svg  # no data points plotted, only legend swatches


@pytest.mark.parametrize("bad_level", [0.0, 1.0, -0.1, 1.5])
def test_rejects_invalid_significance_level(bad_level):
    with pytest.raises(ValueError):
        render_manhattan_svg([], significance_level=bad_level)


def test_real_committed_cohort_scan_run_places_etv6_above_and_ret_braf_below_threshold():
    """Locked to the locus-validated rerun's real numbers. RET's q-value moved
    from ~0.0426 (significant) pre-fix to ~0.1198 (not significant) post-fix
    because the locus validation corrected RET's breakpoint mapping; ETV6 and
    BRAF's significance calls were unaffected by that fix.
    """
    payload = json.loads(REAL_COHORT_SCAN_SUMMARY_JSON.read_text())
    rows = payload["genes"]
    assert len(rows) == 544

    genes_by_symbol = {row["gene_symbol"]: row for row in rows}
    etv6, ret, braf = genes_by_symbol["ETV6"], genes_by_symbol["RET"], genes_by_symbol["BRAF"]
    assert etv6["fdr_significant"] is True
    assert ret["fdr_significant"] is False
    assert braf["fdr_significant"] is False
    assert math.isclose(etv6["min_fdr_adjusted_q_value"], 0.004334, rel_tol=1e-3)
    assert math.isclose(ret["min_fdr_adjusted_q_value"], 0.11976, rel_tol=1e-3)
    assert math.isclose(braf["min_fdr_adjusted_q_value"], 0.23560, rel_tol=1e-3)

    svg = render_manhattan_svg(rows, significance_level=0.05)
    threshold_y = _threshold_line_y(svg)

    etv6_y = _circle_cy(svg, "ETV6")
    ret_y = _circle_cy(svg, "RET")
    braf_y = _circle_cy(svg, "BRAF")

    # SVG y grows downward, so "visually above the dashed line" means a
    # smaller cy than the threshold line's y1/y2.
    assert etv6_y < threshold_y, "ETV6 must plot above the FDR threshold line"
    assert ret_y > threshold_y, "RET must plot below the FDR threshold line"
    assert braf_y > threshold_y, "BRAF must plot below the FDR threshold line"

    # The one significant gene is labeled by symbol.
    assert ">ETV6<" in svg

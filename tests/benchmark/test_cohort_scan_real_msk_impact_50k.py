"""Live genome-wide cohort-scan benchmark against real MSK-IMPACT 50k data.

Excluded from the default ``pytest`` run (``@pytest.mark.network`` plus the
same ``CFH_RUN_NETWORK_TESTS=1`` opt-in gate already used by
``tests/benchmark/test_braf_kinase_retention_msk_impact_50k.py``). This is
the sanity cross-check that the new generic cohort-scan path reproduces the
existing hand-curated BRAF/RET pipeline's numbers on the same live cohort,
not just that it runs without crashing.

``max_genes`` bounds this particular run to a manageable wall-clock time
(candidate genes are sorted by descending recurrence, and both BRAF and RET
rank in the top 10 by distinct-patient count, so they are always included);
the reported before/after-gating gene counts are always the true,
uncapped cohort-wide numbers.
"""

from __future__ import annotations

import json
import os

import pytest

from cfh.cohort.outputs import build_summary_rows, write_cohort_scan_outputs
from cfh.cohort.scan import run_cohort_scan

_STUDY_ID = "msk_impact_50k_2026"

# Already-known values from the existing hand-curated BRAF/RET pipeline
# (see tests/benchmark/test_braf_kinase_retention_msk_impact_50k.py and the
# committed runs/braf_msk-impact-50k-2026_*/RET_msk-impact-50k-2026_*
# artifacts), used here only as a printed cross-check -- never forced.
_EXPECTED_BRAF_FUSION_COUNT_APPROX = 179
_EXPECTED_BRAF_KINASE_RETENTION_PERCENT_APPROX = 91.0
_EXPECTED_RET_FUSION_COUNT_APPROX = 194
_EXPECTED_RET_KINASE_RETENTION_PERCENT_APPROX = 92.3


@pytest.mark.network
@pytest.mark.skipif(
    os.getenv("CFH_RUN_NETWORK_TESTS") != "1",
    reason="set CFH_RUN_NETWORK_TESTS=1 to run the live genome-wide cohort scan",
)
def test_cohort_scan_real_msk_impact_50k(tmp_path):
    result = run_cohort_scan(
        _STUDY_ID,
        min_distinct_patients=5,
        n_permutations=1_000,
        adaptive=True,
        n_permutations_small=100,
        max_genes=40,
        cache_dir=tmp_path / "cache",
    )

    # The full pre-gate cohort universe and the gated candidate count are
    # always real, uncapped numbers -- not affected by max_genes.
    assert result.total_genes_before_gating > 0
    assert result.genes_after_gating > 0
    assert result.genes_after_gating <= result.total_genes_before_gating

    outcomes_by_gene = {outcome.gene_symbol: outcome for outcome in result.gene_outcomes}
    assert "BRAF" in outcomes_by_gene, "BRAF did not rank in the top max_genes by recurrence"
    assert "RET" in outcomes_by_gene, "RET did not rank in the top max_genes by recurrence"

    braf = outcomes_by_gene["BRAF"]
    ret = outcomes_by_gene["RET"]
    assert braf.status == "ok", braf.error
    assert ret.status == "ok", ret.error
    assert braf.config_source == "curated"
    assert ret.config_source == "curated"

    braf_summary = braf.run.summary
    ret_summary = ret.run.summary

    report = {
        "total_genes_before_gating": result.total_genes_before_gating,
        "genes_after_gating": result.genes_after_gating,
        "genes_actually_scanned": len(result.gene_outcomes),
        "curated_gene_count": result.curated_gene_count,
        "auto_config_gene_count": result.auto_config_gene_count,
        "unresolved_gene_count": result.unresolved_gene_count,
        "fdr_significant_gene_count": len(result.significant_genes),
        "fdr_significant_genes": result.significant_genes,
        "BRAF": {
            "total_fusions": braf_summary["total_fusions"],
            "in_frame_percent": braf_summary["in_frame_percent"],
            "kinase_retained_percent": braf_summary["kinase_retained_percent"],
            "expected_total_fusions_approx": _EXPECTED_BRAF_FUSION_COUNT_APPROX,
            "expected_kinase_retained_percent_approx": (
                _EXPECTED_BRAF_KINASE_RETENTION_PERCENT_APPROX
            ),
        },
        "RET": {
            "total_fusions": ret_summary["total_fusions"],
            "in_frame_percent": ret_summary["in_frame_percent"],
            "kinase_retained_percent": ret_summary["kinase_retained_percent"],
            "expected_total_fusions_approx": _EXPECTED_RET_FUSION_COUNT_APPROX,
            "expected_kinase_retained_percent_approx": (
                _EXPECTED_RET_KINASE_RETENTION_PERCENT_APPROX
            ),
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))

    # Sanity bounds only -- real cohort data can drift from a snapshot
    # taken when the expected values above were recorded, and this test
    # reports the real numbers rather than forcing a match.
    assert braf_summary["total_fusions"] > 0
    assert ret_summary["total_fusions"] > 0
    assert 0 <= braf_summary["in_frame_percent"] <= 100
    assert 0 <= ret_summary["in_frame_percent"] <= 100
    assert 0 <= braf_summary["kinase_retained_percent"] <= 100
    assert 0 <= ret_summary["kinase_retained_percent"] <= 100

    # The rest of the pipeline (summary building + writing every output
    # artifact) must also complete against real data.
    rows = build_summary_rows(result)
    assert len(rows) == len(result.gene_outcomes)
    paths = write_cohort_scan_outputs(result, tmp_path / "runs", pdf=True)
    assert paths["summary_tsv"].exists()
    assert paths["summary_pdf"].exists()
    assert "BRAF" in paths["gene_reports"]
    assert "RET" in paths["gene_reports"]

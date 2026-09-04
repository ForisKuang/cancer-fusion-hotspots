"""Live-path benchmark coverage for the original MSK-IMPACT 2017 study."""

from __future__ import annotations

import json
import os

import pytest

from cfh.real_benchmark import run_real_benchmark

_STUDY_ID = "msk_impact_2017"
_PROFILE_ID = "msk_impact_2017_structural_variants"


@pytest.mark.network
@pytest.mark.skipif(
    os.getenv("CFH_RUN_NETWORK_TESTS") != "1",
    reason="set CFH_RUN_NETWORK_TESTS=1 to run live cBioPortal/Genome Nexus benchmark",
)
def test_braf_kinase_retention_in_real_msk_impact_2017():
    """Exercise the shared production benchmark pipeline against the 2017 cohort."""
    run = run_real_benchmark("BRAF", _STUDY_ID, n_permutations=25)

    assert run.study_id == _STUDY_ID
    assert run.molecular_profile_id == _PROFILE_ID
    assert run.raw_structural_variant_count > 0
    assert len(run.events) == len(run.features) > 0
    assert run.summary["mapped_fusions"] > 0
    assert any(event.Frame_status == "in-frame" for event in run.events)
    assert all(
        feature.Domain_retention_flags["kinase"] != "unknown"
        for feature in run.features
    )
    assert 0 <= run.summary["fisher_p_value"] <= 1
    assert sum(
        row["Event_count"] for row in run.summary["partner_counts"]
    ) == run.summary["total_fusions"]
    print(
        json.dumps(
            {
                "study_id": run.study_id,
                "molecular_profile_id": run.molecular_profile_id,
                "raw_structural_variants": run.raw_structural_variant_count,
                "protein_fusions": run.summary["total_fusions"],
                "mapped_fusions": run.summary["mapped_fusions"],
                "in_frame": run.summary["in_frame_count"],
                "kinase_retained": run.summary["kinase_retained_count"],
                "unique_partners": len(run.summary["partner_counts"]),
                "fisher_odds_ratio": run.summary["fisher_odds_ratio"],
                "fisher_p_value": run.summary["fisher_p_value"],
            },
            sort_keys=True,
        )
    )

"""Environment-gated live checks for the WP8 generalization paths."""

from __future__ import annotations

import os

import pytest

from cfh.real_benchmark import run_analysis, run_real_benchmark

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        os.getenv("CFH_RUN_NETWORK_TESTS") != "1",
        reason="set CFH_RUN_NETWORK_TESTS=1 to run live RET/TCGA benchmarks",
    ),
]


def _assert_live_run(run, expected_domain_bounds: tuple[int, int]) -> None:
    assert run.raw_structural_variant_count > 0
    assert len(run.events) == len(run.features) > 0
    assert run.summary["mapped_fusions"] > 0
    assert all(feature.Domain_retention_flags["kinase"] != "unknown" for feature in run.features)
    assert (
        run.summary["domain_start_aa"],
        run.summary["domain_end_aa"],
    ) == expected_domain_bounds


def test_ret_analyze_path_against_real_msk_impact_50k():
    run = run_analysis("RET", "msk_impact_50k_2026", n_permutations=25)

    _assert_live_run(run, (724, 1005))
    assert any(result.Algorithm == "frequency" for result in run.results)


def test_braf_domain_retention_against_real_tcga_pan_cancer_sv():
    run = run_real_benchmark(
        "BRAF",
        "thca_tcga_pan_can_atlas_2018",
        n_permutations=25,
    )

    _assert_live_run(run, (457, 712))
    assert run.endpoints[1].startswith("https://grch38.genomenexus.org/")
    assert all(event.Cohort == "thca_tcga_pan_can_atlas_2018" for event in run.events)

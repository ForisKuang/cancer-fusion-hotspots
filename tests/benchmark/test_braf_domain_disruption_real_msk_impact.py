"""CI-gating BRAF domain-disruption benchmark against already-committed real
MSK-IMPACT run data.

This reuses the checked-in ``runs/braf_msk-impact-*`` directories produced by
the WP10 real-data BRAF benchmark (a live cBioPortal + Genome Nexus pull,
already committed to the repo) instead of re-fetching live data. Each run's
``results.json`` already carries the real per-event breakpoint protein
position and 5'/3' role for every mapped BRAF fusion; this test reclassifies
each event's RAS-binding (PF02196) and cysteine-rich (PF00130) domain status
from those real positions using the production ``classify_domain_retention``
function and the real domain boundaries from the committed Genome Nexus
fixture (the same source used for the kinase domain), then runs the actual
``domain_disruption`` algorithm plugin end to end.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cfh.algorithms.domain_disruption import DomainDisruptionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.mapping.feature_mapper import classify_domain_retention
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RUNS_DIR = _REPO_ROOT / "runs"
_GENOME_NEXUS_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "genome_nexus"
    / "canonical_transcript_braf.json"
)


def _real_run_results_path(prefix: str) -> Path:
    candidates = sorted(_RUNS_DIR.glob(f"{prefix}_*"))
    if not candidates:
        pytest.skip(f"no committed real run directory found for {prefix!r} under {_RUNS_DIR}")
    return candidates[-1] / "results.json"


def _domain_bounds(pfam_id: str) -> tuple[int, int]:
    payload = json.loads(_GENOME_NEXUS_FIXTURE.read_text())
    domain = next(d for d in payload["pfamDomains"] if d["pfamDomainId"] == pfam_id)
    return domain["pfamDomainStart"], domain["pfamDomainEnd"]


def _events_and_features_from_real_run(
    results_path: Path,
) -> tuple[list[FusionEvent], list[FusionFeature]]:
    payload = json.loads(results_path.read_text())
    ras_bounds = _domain_bounds("PF02196")
    cys_bounds = _domain_bounds("PF00130")

    events: list[FusionEvent] = []
    features: list[FusionFeature] = []
    for row in payload["events"]:
        events.append(
            FusionEvent(
                Event_id=row["event_id"],
                Cohort=payload["study_id"],
                Sample_id=row["sample_id"],
                Fusion_name=row["fusion_name"],
                Frame_status=row["frame_status"],
                Is_protein_fusion=True,
            )
        )
        role = row["target_role"]
        position = row["breakpoint_protein_position"]
        features.append(
            FusionFeature(
                Event_id=row["event_id"],
                Gene="BRAF",
                Role=role,
                Junction_position_aa=position,
                Domain_retention_flags={
                    "kinase": row["domain_status"],
                    "ras_binding": classify_domain_retention(*ras_bounds, position, role),
                    "cysteine_rich": classify_domain_retention(*cys_bounds, position, role),
                },
            )
        )
    return events, features


@pytest.mark.parametrize(
    "study_prefix, expected_table, expected_disruption_rate",
    [
        # 2017 cohort (Zehir et al. successor pull): 33 in-frame BRAF fusions.
        ("braf_msk-impact-2017", [[31, 2], [2, 0]], pytest.approx(31 / 33)),
        # 50k 2026 cohort: 151 in-frame BRAF fusions.
        ("braf_msk-impact-50k-2026", [[126, 19], [25, 4]], pytest.approx(126 / 151)),
    ],
)
def test_braf_domain_disruption_against_committed_real_msk_impact_run(
    study_prefix, expected_table, expected_disruption_rate
):
    results_path = _real_run_results_path(study_prefix)
    events, features = _events_and_features_from_real_run(results_path)
    config = load_gene_config("braf")

    result = DomainDisruptionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 5_000}
    )

    print(
        json.dumps(
            {
                "study": study_prefix,
                "n_events": len(events),
                "contingency_table": result.Tables["frame_domain_contingency_table"],
                "observed_in_frame_disruption_rate": result.Summary[
                    "observed_in_frame_disruption_rate"
                ],
                "fisher_odds_ratio": result.Summary["fisher_odds_ratio"],
                "fisher_p_value": result.Summary["fisher_p_value"],
                "permutation_empirical_p_value": result.Summary["permutation_empirical_p_value"],
            },
            sort_keys=True,
        )
    )

    # Real-data reality check: RAS-binding/cysteine-rich domain loss IS the
    # majority pattern among in-frame BRAF fusions (breakpoints cluster past
    # aa ~380, well downstream of both domains at aa 156-280), but it is NOT
    # a clean 100% pattern the way kinase-domain *retention* is -- a handful
    # of real events keep the N-terminal module (BRAF-as-5'-partner
    # readthrough fusions, and a few breakpoints mapped to protein position 1
    # from an intronic clamp). This assertion is intentionally exact rather
    # than a loose bound so a future regression in the classification/stats
    # path is caught, not silently rationalized away.
    assert result.Tables["frame_domain_contingency_table"] == expected_table
    assert result.Summary["observed_in_frame_disruption_rate"] == expected_disruption_rate
    assert 0 <= result.Summary["fisher_p_value"] <= 1
    assert 0 <= result.Summary["permutation_empirical_p_value"] <= 1


def test_domain_disruption_is_not_significantly_enriched_in_frame_in_real_braf_data():
    """Honest negative finding, verified rather than assumed: unlike kinase
    *retention* (which real 2017-cohort data shows IS enriched in-frame,
    p=0.010 per the domain_retention benchmark), N-terminal domain *loss* is
    common across BOTH in-frame and out-of-frame/unknown BRAF fusions in the
    real 2017 cohort -- so it is not a statistically significant in-frame
    enrichment signal by Fisher's exact test, even though the raw rate is
    high (31/33 = 93.9%). Breakpoints cluster well past the domain
    regardless of frame outcome, matching the earlier aa-380+ clustering
    observation; that does not, by itself, imply frame-specific enrichment.
    """
    results_path = _real_run_results_path("braf_msk-impact-2017")
    events, features = _events_and_features_from_real_run(results_path)
    config = load_gene_config("braf")

    result = DomainDisruptionAlgorithm().run(
        events, features, config, {"seed": 42, "n_permutations": 5_000}
    )

    assert result.Summary["observed_in_frame_disruption_rate"] == pytest.approx(31 / 33)
    assert result.Summary["fisher_p_value"] == 1.0
    assert result.Summary["fisher_odds_ratio"] == 0.0

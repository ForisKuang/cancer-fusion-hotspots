"""CI-gating BRAF cutpoint-detection sanity check.

Reuses the domain-retention benchmark's committed fixture and field-mapping
helper (rather than duplicating any cBioPortal ingestion code) to run the
gene-agnostic cutpoint-detection algorithm against a synthetic-but-realistic
fixture -- modeled on the named BRAF-fusion examples in Zehir et al.
(PMC5461196), per that fixture's own source comment -- and checks the
inferred cutpoint against the known kinase-domain boundary (Pfam PF07714, aa
458) parsed from the committed Genome Nexus canonical-transcript fixture.

This is a CI-gating sanity check against synthetic data, not a validation
against real fusion calls: a genuine live-data validation, pulling real
BRAF structural variants through the same cBioPortal + Genome Nexus path the
domain-retention benchmark's `test_braf_kinase_retention_in_real_msk_impact_50k`
(network-marked, opt-in via ``CFH_RUN_NETWORK_TESTS=1``) already exercises,
is a follow-up and has not been done here.
"""

from __future__ import annotations

import json
from pathlib import Path

from test_braf_kinase_retention_msk_impact_50k import _fixture_events_and_features

from cfh.algorithms.cutpoint_detection import CutpointDetectionAlgorithm
from cfh.genes.registry import load_gene_config
from cfh.mapping.genome_nexus_source import parse_canonical_transcript

_GENOME_NEXUS_BRAF_FIXTURE = (
    Path(__file__).parents[1] / "fixtures" / "genome_nexus" / "canonical_transcript_braf.json"
)
_KNOWN_KINASE_DOMAIN_START_AA = 458


def _braf_pfam_domain_boundaries() -> list[int]:
    payload = json.loads(_GENOME_NEXUS_BRAF_FIXTURE.read_text())
    transcript = parse_canonical_transcript(payload)
    boundaries: set[int] = set()
    for domain in transcript.pfam_domains:
        boundaries.add(domain.start_aa)
        boundaries.add(domain.end_aa)
    return sorted(boundaries)


def test_braf_cutpoint_lands_near_known_kinase_domain_boundary():
    events, features = _fixture_events_and_features()
    config = load_gene_config("braf")
    boundaries = _braf_pfam_domain_boundaries()
    assert _KNOWN_KINASE_DOMAIN_START_AA in boundaries

    result = CutpointDetectionAlgorithm().run(
        events,
        features,
        config,
        {"seed": 42, "n_permutations": 1_000, "domain_boundaries": boundaries},
    )

    assert result.Summary["determinable"] is True
    assert result.Summary["corrected_p_value"] < 0.05

    # Sanity check (reported honestly regardless of outcome), against the
    # synthetic fixture described above: its kinase-retained events all sit
    # at or below aa 455 and every kinase-lost/disrupted event sits at aa
    # 500 or above, so the recurrence-only, label-free scan should
    # independently land within a few residues of the known kinase-domain
    # start (aa 458) -- with no BRAF-specific logic anywhere in the
    # algorithm itself.
    comparison = result.Summary["known_domain_boundary_comparison"]
    assert comparison is not None
    assert comparison["nearest_known_domain_boundary_aa"] == _KNOWN_KINASE_DOMAIN_START_AA
    assert comparison["distance_aa"] <= 10

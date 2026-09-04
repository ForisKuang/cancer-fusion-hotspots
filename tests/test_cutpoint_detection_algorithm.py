from unittest.mock import MagicMock

from cfh.algorithms import registry
from cfh.algorithms.cutpoint_detection import (
    CutpointDetectionAlgorithm,
    _known_domain_boundaries,
    _nearest_boundary_comparison,
)
from cfh.genes.registry import GeneConfig, KeyDomain, load_gene_config
from cfh.mapping.domain_source import ProteinDomain
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

# A fabricated single-gene config, deliberately unrelated to any real gene
# in genes/configs/, to demonstrate the algorithm has no gene-specific
# hardcoding of its own.
_FAKE_GENE = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    key_domains=[KeyDomain(name="Made-up domain", source="test", key="made_up")],
)


def _events_and_features(records: list[tuple[str, int, str]]) -> tuple[list, list]:
    events = []
    features = []
    for event_id, position, status in records:
        events.append(FusionEvent(Event_id=event_id, Cohort="synthetic"))
        features.append(
            FusionFeature(
                Event_id=event_id,
                Gene="FAKE1",
                Junction_position_aa=position,
                Domain_retention_flags={"made_up": status},
            )
        )
    return events, features


def _clean_separation_records() -> list[tuple[str, int, str]]:
    retained = [("r0", 100, "retained"), ("r1", 110, "retained"), ("r2", 120, "retained")]
    lost = [("l0", 300, "lost"), ("l1", 310, "disrupted"), ("l2", 320, "lost")]
    return retained + lost


def test_algorithm_registered():
    assert registry.get("cutpoint_detection") is CutpointDetectionAlgorithm
    assert "cutpoint_detection" in registry.list_algorithms()


def test_recovers_cutpoint_from_synthetic_clean_separation():
    events, features = _events_and_features(_clean_separation_records())

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 300}
    )

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "cutpoint_detection"
    assert result.Summary["determinable"] is True
    assert result.Summary["inferred_cutpoint_aa"] == 120
    assert result.Summary["n_events_analyzed"] == 6
    assert 0.0 <= result.Summary["corrected_p_value"] <= 1.0
    assert result.Warnings == []
    assert result.Tables["cutpoint_scan"]
    assert result.Summary["known_domain_boundary_comparison"] is None


def test_domain_boundary_comparison_reports_nearest_distance():
    events, features = _events_and_features(_clean_separation_records())

    result = CutpointDetectionAlgorithm().run(
        events,
        features,
        _FAKE_GENE,
        {"seed": 42, "n_permutations": 300, "domain_boundaries": [118, 400]},
    )

    comparison = result.Summary["known_domain_boundary_comparison"]
    assert comparison == {"nearest_known_domain_boundary_aa": 118, "distance_aa": 2}


def test_domain_boundary_distance_is_signed_and_reported_outside_every_domain_span():
    """The inferred cutpoint falling outside every configured domain span
    (issue #14) must still report the nearest boundary and a *signed*
    distance, not ``None``: positive when the cutpoint is beyond (C-terminal
    of) the nearest boundary, negative when it's before (N-terminal of) it.
    """
    assert _nearest_boundary_comparison(120, [50, 90]) == {
        "nearest_known_domain_boundary_aa": 90,
        "distance_aa": 30,
    }
    assert _nearest_boundary_comparison(10, [50, 90]) == {
        "nearest_known_domain_boundary_aa": 50,
        "distance_aa": -40,
    }


def test_ret_real_committed_cutpoint_gets_nearest_boundary_not_null_without_gn_pfam_domains():
    """Regression test for issue #14.

    A real full-cohort run for RET (``runs/ret_msk-impact-50k-2026_20260903T160639Z/
    results.json``) inferred a cutpoint of 1063 aa from 194 real events and
    recorded ``known_domain_boundary_comparison: null`` -- Genome Nexus's
    canonical-transcript response for RET carried no Pfam domain annotation
    at the time, so the old genome-nexus-only lookup came back empty and the
    comparison was silently dropped, even though 1063 aa is a well-defined
    58 aa past the real kinase domain's end (1005 aa, PF07714, verified live
    against https://www.genomenexus.org and UniProt P07949) -- i.e. outside
    every known domain span, not "no domains known at all". Every other
    committed real run for RET (e.g. the canonical full-cohort run at
    ``runs/ret_msk-impact-50k-2026_20260903T193620Z/results.json``, and
    ``runs/cohort-scan_msk_impact_50k_2026_20260903T223024Z/cohort_scan/
    gene_reports/ret/results.json``) recovers a boundary comparison for this
    same 1063 aa/1005 aa pair, but with an *unsigned* ``distance_aa`` (``58``
    either way it's approached) -- this test also proves the distance is now
    genuinely signed (see
    ``test_domain_boundary_distance_is_signed_and_reported_outside_every_domain_span``
    for the negative-sign case).

    This reproduces the historical Genome-Nexus-has-nothing case explicitly
    (an empty ``pfamDomains`` payload) and proves the UniProt fallback (the
    same one ``domain_retention`` already relies on via
    :func:`~cfh.mapping.genome_nexus_source.resolve_domains`) now recovers
    that real, correct, non-null comparison instead.
    """
    ret_config = load_gene_config("RET")
    real_committed_cutpoint_aa = 1063

    empty_genome_nexus_client = MagicMock()
    empty_genome_nexus_client.fetch_canonical_transcript.return_value = {
        "transcriptId": ret_config.canonical_transcript_id,
        "pfamDomains": [],
    }
    real_ret_uniprot_domains = MagicMock()
    real_ret_uniprot_domains.fetch.return_value = [
        ProteinDomain(name="PF00028", start_aa=173, end_aa=256, source="uniprot"),
        ProteinDomain(name="Protein kinase domain", start_aa=724, end_aa=1005, source="uniprot"),
    ]

    boundaries = _known_domain_boundaries(
        ret_config,
        {
            "genome_nexus_client": empty_genome_nexus_client,
            "uniprot_source": real_ret_uniprot_domains,
        },
    )

    assert boundaries == [173, 256, 724, 1005]
    assert _nearest_boundary_comparison(real_committed_cutpoint_aa, boundaries) == {
        "nearest_known_domain_boundary_aa": 1005,
        "distance_aa": 58,
    }
    real_ret_uniprot_domains.fetch.assert_called_once_with(ret_config.protein_id)


def test_degenerate_input_is_reported_as_not_determinable_without_raising():
    events, features = _events_and_features(
        [("a", 100, "retained"), ("b", 110, "retained"), ("c", 120, "retained")]
    )

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["determinable"] is False
    assert result.Summary["inferred_cutpoint_aa"] is None
    assert result.Warnings
    assert result.Summary["known_domain_boundary_comparison"] is None


def test_events_missing_mapped_positions_are_excluded_not_crashed_on():
    events = [FusionEvent(Event_id="unmapped", Cohort="synthetic")]
    features = [
        FusionFeature(
            Event_id="unmapped",
            Gene="FAKE1",
            Junction_position_aa=None,
            Domain_retention_flags={"made_up": "retained"},
        )
    ]

    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )

    assert result.Summary["n_events_analyzed"] == 0
    assert result.Summary["determinable"] is False


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _events_and_features(_clean_separation_records())
    result = CutpointDetectionAlgorithm().run(
        events, features, _FAKE_GENE, {"seed": 42, "n_permutations": 50}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)

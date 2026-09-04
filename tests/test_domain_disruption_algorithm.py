import pytest

from cfh.algorithms import registry
from cfh.algorithms.domain_disruption import DomainDisruptionAlgorithm
from cfh.genes.registry import GeneConfig, KeyDomain, load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import DomainRetentionDetail, FusionFeature
from cfh.stats.breakpoint_tests import (
    DISRUPTED_STATUSES,
    build_frame_domain_contingency_table,
    fishers_frame_domain_test,
    permutation_null_test,
)

# A fabricated single-gene config, deliberately unrelated to any real gene
# in genes/configs/, to demonstrate the algorithm has no gene-specific
# hardcoding of its own.
_FAKE_GENE_WITH_DISRUPTION_DOMAIN = GeneConfig(
    gene_symbol="FAKE1",
    canonical_transcript_id="NM_000001",
    protein_id="P00001",
    disruption_required_domains=[
        KeyDomain(name="Made-up autoinhibitory domain", source="test", key="made_up")
    ],
)

_FAKE_GENE_NO_DISRUPTION_DOMAIN = GeneConfig(
    gene_symbol="FAKE2",
    canonical_transcript_id="NM_000002",
    protein_id="P00002",
    key_domains=[KeyDomain(name="Some other domain", source="test", key="other")],
)


def _events_and_features(records: list[tuple[str, str, str, int]]) -> tuple[list, list]:
    """Build (event, feature) pairs from ``(event_id, frame_status, domain_status,
    breakpoint_position)``."""
    events = []
    features = []
    for event_id, frame_status, status, position in records:
        events.append(
            FusionEvent(
                Event_id=event_id,
                Cohort="synthetic",
                Frame_status=frame_status,
                Is_protein_fusion=True,
            )
        )
        features.append(
            FusionFeature(
                Event_id=event_id,
                Gene="FAKE1",
                Junction_position_aa=position,
                Domain_retention_flags={"made_up": status},
            )
        )
    return events, features


def _clean_separation_records() -> list[tuple[str, str, str, int]]:
    """In-frame fusions all lose the domain; other fusions all retain it."""
    disrupted_in_frame = [
        ("d0", "in-frame", "lost", 300),
        ("d1", "in-frame", "disrupted", 310),
        ("d2", "in-frame", "lost", 320),
        ("d3", "in-frame", "lost", 330),
        ("d4", "in-frame", "disrupted", 340),
        ("d5", "in-frame", "lost", 350),
    ]
    retained_other = [
        ("r0", "out-of-frame", "retained", 100),
        ("r1", "out-of-frame", "retained", 110),
        ("r2", "unknown", "retained", 120),
        ("r3", "unknown", "retained", 130),
    ]
    return disrupted_in_frame + retained_other


def test_algorithm_registered():
    assert registry.get("domain_disruption") is DomainDisruptionAlgorithm
    assert "domain_disruption" in registry.list_algorithms()


def test_result_schema_matches_canonical_algorithm_result_fields():
    events, features = _events_and_features(_clean_separation_records())
    result = DomainDisruptionAlgorithm().run(
        events, features, _FAKE_GENE_WITH_DISRUPTION_DOMAIN, {"seed": 42, "n_permutations": 200}
    )
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_detects_enriched_disruption_in_synthetic_clean_separation():
    events, features = _events_and_features(_clean_separation_records())

    table = build_frame_domain_contingency_table(
        events,
        features,
        _FAKE_GENE_WITH_DISRUPTION_DOMAIN,
        domains=_FAKE_GENE_WITH_DISRUPTION_DOMAIN.disruption_required_domains,
        hit_statuses=DISRUPTED_STATUSES,
    )
    assert table == [[6, 0], [0, 4]]
    odds_ratio, p_value = fishers_frame_domain_test(table)
    assert odds_ratio == float("inf")

    result = DomainDisruptionAlgorithm().run(
        events, features, _FAKE_GENE_WITH_DISRUPTION_DOMAIN, {"seed": 42, "n_permutations": 2_000}
    )

    assert result.Algorithm == "domain_disruption"
    assert result.Tables["frame_domain_contingency_table"] == table
    assert result.Summary["fisher_odds_ratio"] == odds_ratio
    assert result.Summary["fisher_p_value"] == p_value
    assert result.Summary["fisher_p_value"] < 0.05
    assert result.Summary["observed_in_frame_disruption_rate"] == 1.0


def test_disruption_result_reports_optional_truncation_descriptives():
    events, features = _events_and_features(_clean_separation_records())
    features[0] = features[0].model_copy(
        update={
            "Domain_retention_details": {
                "made_up": DomainRetentionDetail(
                    Domain_start_aa=100,
                    Domain_end_aa=200,
                    Retained_start_aa=100,
                    Retained_end_aa=160,
                    Retained_fraction=61 / 101,
                    Is_truncated=True,
                )
            }
        }
    )
    features[1] = features[1].model_copy(
        update={
            "Domain_retention_details": {
                "made_up": DomainRetentionDetail(
                    Domain_start_aa=100,
                    Domain_end_aa=200,
                    Retained_fraction=0.0,
                    Is_truncated=False,
                )
            }
        }
    )

    result = DomainDisruptionAlgorithm().run(
        events,
        features,
        _FAKE_GENE_WITH_DISRUPTION_DOMAIN,
        {"seed": 42, "n_permutations": 10},
    )
    row = result.Tables["domain_disruption_descriptives"][0]

    assert row["Quantitative_call_count"] == 2
    assert row["Truncated_count"] == 1
    assert row["Fully_lost_count"] == 1
    assert row["Mean_retained_fraction_among_non_retained_calls"] == pytest.approx((61 / 101) / 2)


def test_permutation_null_is_bit_identical_for_a_fixed_seed():
    events, features = _events_and_features(_clean_separation_records())
    domains = _FAKE_GENE_WITH_DISRUPTION_DOMAIN.disruption_required_domains
    first = permutation_null_test(
        events,
        features,
        _FAKE_GENE_WITH_DISRUPTION_DOMAIN,
        seed=1729,
        n_permutations=250,
        domains=domains,
        hit_statuses=DISRUPTED_STATUSES,
    )
    second = permutation_null_test(
        events,
        features,
        _FAKE_GENE_WITH_DISRUPTION_DOMAIN,
        seed=1729,
        n_permutations=250,
        domains=domains,
        hit_statuses=DISRUPTED_STATUSES,
    )
    assert first == second


def test_no_op_when_gene_leaves_disruption_required_domains_unset():
    """Same graceful-skip pattern already proven for exon_retention/joint_partner/
    confidence_stats: a gene that never opts in produces a documented no-op,
    not an error."""
    result = DomainDisruptionAlgorithm().run([], [], _FAKE_GENE_NO_DISRUPTION_DOMAIN, {})

    assert result.Algorithm == "domain_disruption"
    assert result.Summary["fisher_p_value"] is None
    assert result.Summary["fisher_odds_ratio"] is None
    assert result.Summary["permutation_empirical_p_value"] is None
    assert result.Summary["observed_in_frame_disruption_rate"] is None
    assert result.Tables == {}
    assert result.Warnings
    assert "FAKE2" in result.Warnings[0]
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_no_op_for_eml4_alk_negative_control():
    """EML4-ALK's oncogenicity is not modeled as N-terminal domain excision; the
    real config deliberately leaves disruption_required_domains unset."""
    config = load_gene_config("eml4-alk")
    assert config.disruption_required_domains == []

    result = DomainDisruptionAlgorithm().run([], [], config, {})

    assert result.Summary["fisher_p_value"] is None
    assert result.Warnings


def test_no_op_for_tmprss2_erg_negative_control():
    """TMPRSS2-ERG's mechanism is promoter-swap/expression-driven, not N-terminal
    domain disruption -- this must no-op rather than report a spurious
    disruption-requirement signal for a gene pair with no such config."""
    config = GeneConfig(gene_pair=("TMPRSS2", "ERG"), analysis_modes=["promoter_swap"])
    assert config.disruption_required_domains == []

    events = [
        FusionEvent(
            Event_id="tmprss2-erg-1",
            Cohort="synthetic",
            Frame_status="in-frame",
            Is_protein_fusion=True,
        )
    ]
    features = [
        FusionFeature(
            Event_id="tmprss2-erg-1",
            Gene="ERG",
            Domain_retention_flags={"ets_domain": "retained"},
        )
    ]

    result = DomainDisruptionAlgorithm().run(events, features, config, {})

    assert result.Summary["fisher_p_value"] is None
    assert result.Summary["observed_in_frame_disruption_rate"] is None
    assert any("disruption_required_domains" in warning for warning in result.Warnings)

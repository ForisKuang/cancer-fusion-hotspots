import pytest

from cfh.algorithms import registry
from cfh.algorithms.joint_partner import JointPartnerAlgorithm, JointPartnerMode
from cfh.genes.registry import GeneConfig, load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent


def _build_synthetic_cohort() -> list[FusionEvent]:
    """Construct a cohort of 100 events where:
    - (GENE_A, GENE_B) is heavily enriched over independence:
      marginal A = 20, marginal B = 20 -> independence expected = 4.0, observed = 18.
    - (GENE_C, GENE_D) occurs at independence rate:
      marginal C = 20, marginal D = 20 -> independence expected = 4.0, observed = 4.
    - Other background events make up the rest of the cohort.
    """
    events: list[FusionEvent] = []
    idx = 0

    # Overrepresented pair: (GENE_A, GENE_B) observed 18 times
    for _ in range(18):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene="GENE_A",
                Three_prime_gene="GENE_B",
            )
        )
        idx += 1

    # Remaining marginal 5' GENE_A pairings with other 3' genes (2 times)
    for i in range(2):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene="GENE_A",
                Three_prime_gene=f"BG_3P_{i}",
            )
        )
        idx += 1

    # Remaining marginal 3' GENE_B pairings with other 5' genes (2 times)
    for i in range(2):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene=f"BG_5P_{i}",
                Three_prime_gene="GENE_B",
            )
        )
        idx += 1

    # Independence-rate pair: (GENE_C, GENE_D) observed 4 times
    for _ in range(4):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene="GENE_C",
                Three_prime_gene="GENE_D",
            )
        )
        idx += 1

    # Remaining marginal 5' GENE_C pairings (16 times)
    for i in range(16):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene="GENE_C",
                Three_prime_gene=f"BG_3P_OTHER_{i}",
            )
        )
        idx += 1

    # Remaining marginal 3' GENE_D pairings (16 times)
    for i in range(16):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene=f"BG_5P_OTHER_{i}",
                Three_prime_gene="GENE_D",
            )
        )
        idx += 1

    # Remaining unlinked background pairs to reach 100 events total (100 - 58 = 42)
    for i in range(42):
        events.append(
            FusionEvent(
                Event_id=f"ev_{idx}",
                Cohort="synthetic",
                Five_prime_gene=f"BG_5P_NULL_{i}",
                Three_prime_gene=f"BG_3P_NULL_{i}",
            )
        )
        idx += 1

    assert len(events) == 100
    return events


def test_algorithm_registration():
    """Verify algorithm is registered via the standard registry mechanism."""
    cls = registry.get("joint_partner")
    assert cls is JointPartnerAlgorithm
    assert "joint_partner" in registry.list_algorithms()


def test_positive_and_negative_control_enrichment():
    """Test positive (overrepresented) and negative (null-rate) pairs in ONE test.

    Contract:
    - Build one synthetic fixture event list where a configured pair (gene_a, gene_b)
      is deliberately overrepresented relative to what each gene's marginal partner
      distribution would predict under independence, AND a second pair (gene_c, gene_d)
      that occurs at close to the independence-predicted rate.
    - In ONE test: assert the configured/overrepresented pair's result has p_value < 0.05,
      and assert the other pair is either absent from results or has p_value >= 0.05
      (both assertions, same test, positive+negative control).
    """
    events = _build_synthetic_cohort()
    algo = JointPartnerAlgorithm()

    # Positive control: configured overrepresented pair (GENE_A, GENE_B)
    config_ab = GeneConfig(
        gene_symbol="GENE_A-GENE_B",
        gene_pair=["GENE_A", "GENE_B"],
    )
    result_ab = algo.run(events=events, features=[], gene_config=config_ab, params={})

    # Assert overrepresented pair has p_value < 0.05
    assert result_ab.Summary is not None
    assert result_ab.Summary["p_value"] < 0.05
    assert result_ab.Summary["is_significant"] is True
    assert result_ab.Summary["observed_count"] == 18
    assert result_ab.Summary["expected_count"] == pytest.approx(4.0)

    # Negative control: pair (GENE_C, GENE_D) at independence rate
    # Checked via Tables of evaluated pairs in the same result:
    assert result_ab.Tables is not None
    pair_stats = result_ab.Tables["pair_stats"]
    cd_row = next(
        (r for r in pair_stats if r["gene_5p"] == "GENE_C" and r["gene_3p"] == "GENE_D"),
        None,
    )
    # Either absent or has p_value >= 0.05:
    if cd_row is not None:
        assert cd_row["p_value"] >= 0.05
        assert cd_row["is_significant"] is False

    # Also verify when (GENE_C, GENE_D) is explicitly configured as target:
    config_cd = GeneConfig(
        gene_symbol="GENE_C-GENE_D",
        gene_pair=["GENE_C", "GENE_D"],
    )
    result_cd = algo.run(events=events, features=[], gene_config=config_cd, params={})
    assert result_cd.Summary is not None
    assert result_cd.Summary["p_value"] >= 0.05
    assert result_cd.Summary["is_significant"] is False

    # Also verify an unobserved / absent pair is absent from pair_stats
    absent_row = next(
        (r for r in pair_stats if r["gene_5p"] == "NONEXISTENT_A"),
        None,
    )
    assert absent_row is None


def test_runs_with_empty_features_and_no_domain_info():
    """Test JointPartnerMode().run succeeds with features=[] and no domain data."""
    events = [
        FusionEvent(
            Event_id="e1",
            Cohort="c1",
            Five_prime_gene="EML4",
            Three_prime_gene="ALK",
        ),
        FusionEvent(
            Event_id="e2",
            Cohort="c1",
            Five_prime_gene="EML4",
            Three_prime_gene="ALK",
        ),
    ]

    # Load EML4-ALK YAML config (has no transcript or domain boundaries)
    config = load_gene_config("EML4-ALK")
    assert config.key_domains == []

    mode = JointPartnerMode()
    # Explicitly pass features=[] and gene_config with no domain info
    result = mode.run(events=events, features=[], gene_config=config, params={})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "joint_partner"
    assert result.Summary is not None
    assert result.Summary["observed_count"] == 2


def test_algorithm_result_schema_matches_canonical():
    """Assert exact field set of AlgorithmResult matches canonical model_fields."""
    events = [
        FusionEvent(
            Event_id="e1",
            Cohort="c1",
            Five_prime_gene="EML4",
            Three_prime_gene="ALK",
        )
    ]
    algo = JointPartnerAlgorithm()
    config = load_gene_config("EML4-ALK")
    result = algo.run(events=events, features=[], gene_config=config, params={})

    # Exact field set comparison with canonical AlgorithmResult
    assert set(result.model_fields.keys()) == set(AlgorithmResult.model_fields.keys())


def test_runs_with_none_config_and_empty_events():
    """Test robust execution when gene_config is None and/or events is empty."""
    algo = JointPartnerAlgorithm()
    result = algo.run(events=[], features=[], gene_config=None, params={})
    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "joint_partner"
    assert result.Summary["total_events"] == 0
    assert result.Warnings == ["No fusion events provided to joint_partner analysis"]
    assert set(result.model_fields.keys()) == set(AlgorithmResult.model_fields.keys())

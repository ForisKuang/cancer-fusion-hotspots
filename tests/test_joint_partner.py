from cfh.algorithms import registry
from cfh.algorithms.joint_partner import JointPartnerMode
from cfh.genes.registry import GeneConfig, load_gene_config
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.stats.joint_partner_stats import calculate_pair_enrichment


def _event(event_id: int, gene5: str, gene3: str) -> FusionEvent:
    return FusionEvent(
        Event_id=f"event-{event_id}",
        Cohort="synthetic",
        Five_prime_gene=gene5,
        Three_prime_gene=gene3,
    )


def _repeated_events(start: int, count: int, gene5: str, gene3: str) -> list[FusionEvent]:
    return [_event(start + index, gene5, gene3) for index in range(count)]


def test_joint_partner_detects_enriched_pair_but_not_independence_control():
    events = [
        *_repeated_events(0, 50, "EML4", "ALK"),
        *_repeated_events(50, 50, "EML4", "OTHER_5"),
        *_repeated_events(100, 50, "OTHER_3", "ALK"),
        *_repeated_events(150, 10, "GENE_C", "GENE_D"),
        *_repeated_events(160, 90, "GENE_C", "OTHER_C"),
        *_repeated_events(250, 90, "OTHER_D", "GENE_D"),
        *_repeated_events(340, 660, "BACKGROUND_5", "BACKGROUND_3"),
    ]
    config = load_gene_config("eml4-alk")

    result = JointPartnerMode().run(events, features=[], gene_config=config, params={})
    configured_pair = result.Tables["pair_results"][0]
    control_pair = calculate_pair_enrichment(events, "GENE_C", "GENE_D")

    assert configured_pair["p_value"] < 0.05
    assert control_pair.p_value >= 0.05


def test_joint_partner_uses_common_plugin_registry_and_needs_no_domain_data():
    config = GeneConfig(gene_pair=("EML4", "ALK"), analysis_modes=["joint_partner_dependency"])
    result = JointPartnerMode().run(
        [
            FusionEvent(
                Event_id="site-order-event",
                Cohort="synthetic",
                Site1_gene="EML4",
                Site2_gene="ALK",
            )
        ],
        features=[],
        gene_config=config,
        params={},
    )

    assert registry.get("joint_partner") is JointPartnerMode
    assert result.Tables["pair_results"][0]["observed_count"] == 1
    assert set(type(result).model_fields) == set(AlgorithmResult.model_fields)


def test_pair_config_loads_through_gene_registry():
    config = load_gene_config("eml4-alk")

    assert config.gene_pair == ("EML4", "ALK")
    assert config.key_domains == []


def test_gene_config_with_no_gene_pair_is_a_clean_noop_not_a_raise():
    """A gene config lacking the required ``gene_pair`` field (e.g. a real
    full-cohort run for a single-gene config like RET/BRAF that never opted
    into this optional pair-enrichment analysis) must produce a no-op
    result, not an ``Algorithm failed`` warning from a raised exception.
    """
    config = load_gene_config("RET")

    result = JointPartnerMode().run([], features=[], gene_config=config, params={})

    assert isinstance(result, AlgorithmResult)
    assert result.Algorithm == "joint_partner"
    assert result.Summary == {}
    assert result.Warnings == [
        "RET has no gene_pair configured; joint-partner analysis was skipped."
    ]

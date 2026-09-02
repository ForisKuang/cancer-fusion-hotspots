"""Frame/domain-retention benchmark algorithm plugin."""

from __future__ import annotations

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import (
    build_frame_domain_contingency_table,
    fishers_frame_domain_test,
    permutation_null_test,
)


@register("domain_retention")
class DomainRetentionAlgorithm(Algorithm):
    """Test whether configured-domain retention is enriched in-frame."""

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        table = build_frame_domain_contingency_table(events, features, gene_config)
        odds_ratio, fisher_p_value = fishers_frame_domain_test(table)
        permutation_p_value, observed_rate, null_rates = permutation_null_test(
            events,
            features,
            gene_config,
            seed=params.get("seed", 42),
            n_permutations=params.get("n_permutations", 10_000),
            genome_nexus_client=params.get("genome_nexus_client"),
        )
        return AlgorithmResult(
            Algorithm="domain_retention",
            Algorithm_version="0.1.0",
            Parameters={
                "seed": params.get("seed", 42),
                "n_permutations": params.get("n_permutations", 10_000),
            },
            Summary={
                "fisher_odds_ratio": odds_ratio,
                "fisher_p_value": fisher_p_value,
                "permutation_empirical_p_value": permutation_p_value,
                "observed_in_frame_retention_rate": observed_rate,
            },
            Tables={
                "frame_domain_contingency_table": table,
                "permutation_null_retention_rates": list(null_rates),
            },
        )

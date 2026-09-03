"""Frame/domain-disruption benchmark algorithm plugin -- the inverse of
:mod:`cfh.algorithms.domain_retention`.

Some oncogenic fusions require a configured domain to be LOST/DISRUPTED
rather than retained (e.g. an autoinhibitory regulatory module excised by
the breakpoint). This reuses the exact same Fisher's-exact/permutation
machinery as domain retention, just testing domain *exclusion* instead of
domain *retention* as the enriched outcome among in-frame fusions.

Opt-in per gene via ``GeneConfig.disruption_required_domains``: a gene that
leaves this field unset (its default, empty list) produces a no-op result
with no statistics computed, the same graceful-skip pattern already used by
``exon_retention``/``joint_partner``/``confidence_stats`` for genes that
don't configure their respective optional fields.
"""

from __future__ import annotations

from datetime import datetime, timezone

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import (
    DISRUPTED_STATUSES,
    build_frame_domain_contingency_table,
    domain_retention_descriptive_table,
    fishers_frame_domain_test,
    permutation_null_test,
)

VERSION = "0.2.0"


@register("domain_disruption")
class DomainDisruptionAlgorithm(Algorithm):
    """Test whether configured-domain loss/disruption is enriched in-frame."""

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        domains = gene_config.disruption_required_domains
        if not domains:
            return AlgorithmResult(
                Algorithm="domain_disruption",
                Algorithm_version=VERSION,
                Parameters={},
                Summary={
                    "fisher_odds_ratio": None,
                    "fisher_p_value": None,
                    "permutation_empirical_p_value": None,
                    "observed_in_frame_disruption_rate": None,
                },
                Tables={},
                Warnings=[
                    f"{gene_config.gene_symbol or gene_config.gene_pair} has no "
                    "disruption_required_domains configured; domain-disruption "
                    "analysis was skipped."
                ],
                Created_at=datetime.now(timezone.utc),
            )

        table = build_frame_domain_contingency_table(
            events, features, gene_config, domains=domains, hit_statuses=DISRUPTED_STATUSES
        )
        odds_ratio, fisher_p_value = fishers_frame_domain_test(table)
        permutation_p_value, observed_rate, null_rates = permutation_null_test(
            events,
            features,
            gene_config,
            domains=domains,
            hit_statuses=DISRUPTED_STATUSES,
            seed=params.get("seed", 42),
            n_permutations=params.get("n_permutations", 10_000),
            genome_nexus_client=params.get("genome_nexus_client"),
        )
        return AlgorithmResult(
            Algorithm="domain_disruption",
            Algorithm_version=VERSION,
            Parameters={
                "seed": params.get("seed", 42),
                "n_permutations": params.get("n_permutations", 10_000),
            },
            Summary={
                "fisher_odds_ratio": odds_ratio,
                "fisher_p_value": fisher_p_value,
                "permutation_empirical_p_value": permutation_p_value,
                "observed_in_frame_disruption_rate": observed_rate,
            },
            Tables={
                "frame_domain_contingency_table": table,
                "domain_disruption_descriptives": domain_retention_descriptive_table(
                    features, gene_config, domains=domains
                ),
                "permutation_null_disruption_rates": list(null_rates),
            },
        )

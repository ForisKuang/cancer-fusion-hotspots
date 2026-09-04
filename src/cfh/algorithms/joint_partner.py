"""Joint-partner-dependent fusion-oncogenicity analysis mode.

Opt-in via ``GeneConfig.gene_pair``: a gene config with no ``gene_pair`` set
produces a no-op result with no statistics computed, the same graceful-skip
pattern already used by ``exon_retention``/``domain_disruption`` for genes
that don't configure their respective optional fields.
"""

from __future__ import annotations

from datetime import datetime, timezone

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.joint_partner_stats import calculate_pair_enrichment


@register("joint_partner")
class JointPartnerMode(Algorithm):
    """Assess whether a configured 5'/3' gene pair is jointly enriched.

    This mode deliberately uses only gene-pair identities from ``FusionEvent``;
    domain-level ``FusionFeature`` data is not part of its first-pass model.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        """Run pair enrichment for ``gene_config.gene_pair``."""
        del features
        if gene_config.gene_pair is None:
            return AlgorithmResult(
                Algorithm="joint_partner",
                Algorithm_version="0.1.0",
                Parameters={},
                Summary={},
                Tables={},
                Warnings=[
                    f"{gene_config.gene_symbol or 'This gene'} has no gene_pair "
                    "configured; joint-partner analysis was skipped."
                ],
                Created_at=datetime.now(timezone.utc),
            )

        gene5, gene3 = gene_config.gene_pair
        enrichment = calculate_pair_enrichment(events, gene5, gene3)
        significance_level = float(params.get("significance_level", 0.05))
        warnings = []
        if enrichment.eligible_event_count == 0:
            warnings.append("No events with a complete gene pair were available for testing.")

        return AlgorithmResult(
            Algorithm="joint_partner",
            Algorithm_version="0.1.0",
            Parameters={"significance_level": significance_level},
            Summary={
                "configured_pair": [gene5, gene3],
                "observed_count": enrichment.observed_count,
                "expected_count": enrichment.expected_count,
                "p_value": enrichment.p_value,
                "is_enriched": enrichment.p_value < significance_level,
            },
            Tables={"pair_results": [enrichment.as_dict()]},
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )

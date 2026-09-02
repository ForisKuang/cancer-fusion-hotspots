"""Joint-partner-dependent fusion hotspot detection algorithm.

Unlike single-gene domain-retention algorithms, this algorithm evaluates
oncogenicity driven by the specific combination of two partner genes
(e.g., EML4-ALK where EML4's coiled-coil domain enables constitutive
dimerization of ALK's kinase domain).

This module requires only FusionEvent partner annotations (Site1_gene/Site2_gene
or Five_prime_gene/Three_prime_gene) and operates cleanly without FusionFeature
or domain data (i.e. features=[]).
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.joint_partner_stats import evaluate_all_pairs, evaluate_gene_pair

ALGORITHM_NAME = "joint_partner"
ALGORITHM_VERSION = "1.0.0"


def _resolve_target_pair(
    gene_config: Optional[GeneConfig],
    params: dict[str, Any],
) -> Optional[tuple[str, str]]:
    """Determine the target (5', 3') gene pair from params or gene_config."""
    # 1. params overrides
    if "gene_pair" in params and params["gene_pair"]:
        pair = params["gene_pair"]
        if len(pair) >= 2:
            return (str(pair[0]).strip().upper(), str(pair[1]).strip().upper())
    if params.get("gene_5p") and params.get("gene_3p"):
        return (str(params["gene_5p"]).strip().upper(), str(params["gene_3p"]).strip().upper())

    # 2. gene_config
    if gene_config is not None:
        if gene_config.gene_pair and len(gene_config.gene_pair) >= 2:
            return (
                str(gene_config.gene_pair[0]).strip().upper(),
                str(gene_config.gene_pair[1]).strip().upper(),
            )
        if gene_config.partner_5p and gene_config.partner_3p:
            return (
                str(gene_config.partner_5p).strip().upper(),
                str(gene_config.partner_3p).strip().upper(),
            )
        if "-" in gene_config.gene_symbol:
            parts = gene_config.gene_symbol.split("-", 1)
            return (parts[0].strip().upper(), parts[1].strip().upper())

    return None


@register(ALGORITHM_NAME)
@register("joint_partner_dependency")
class JointPartnerAlgorithm(Algorithm):
    """Joint-partner-dependent fusion hotspot detection algorithm.

    Tests whether observed gene pairs occur significantly more frequently
    than expected under an independence null model based on marginal partner
    frequencies in the cohort.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: Optional[GeneConfig],
        params: Optional[dict[str, Any]] = None,
    ) -> AlgorithmResult:
        """Run joint-partner co-occurrence analysis on fusion events.

        Args:
            events: List of normalized FusionEvent objects.
            features: List of FusionFeature objects (unused by joint-partner mode,
                safely accepts empty list).
            gene_config: Optional GeneConfig (e.g. for EML4-ALK). Does not require
                any domain or transcript annotations.
            params: Optional algorithm parameters:
                - directional (bool, default True): whether gene order matters
                - method (str, default 'fisher'): 'fisher' or 'binomial'
                - alpha (float, default 0.05): significance threshold
                - gene_pair (tuple/list of 2 str, optional): target pair override

        Returns:
            AlgorithmResult containing pair co-occurrence statistics.
        """
        params = params or {}
        directional = bool(params.get("directional", True))
        method = str(params.get("method", "fisher"))
        alpha = float(params.get("alpha", 0.05))

        warnings: list[str] = []
        if not events:
            warnings.append("No fusion events provided to joint_partner analysis")

        target_pair = _resolve_target_pair(gene_config, params)

        all_pair_stats = evaluate_all_pairs(
            events,
            directional=directional,
            method=method,
            alpha=alpha,
        )

        summary: dict[str, Any] = {
            "algorithm": ALGORITHM_NAME,
            "version": ALGORITHM_VERSION,
            "directional": directional,
            "method": method,
            "alpha": alpha,
            "total_events": len(events),
            "n_pairs_evaluated": len(all_pair_stats),
        }

        if target_pair is not None:
            target_stats = evaluate_gene_pair(
                events,
                gene_5p=target_pair[0],
                gene_3p=target_pair[1],
                directional=directional,
                method=method,
                alpha=alpha,
            )
            summary.update({
                "gene_5p": target_stats["gene_5p"],
                "gene_3p": target_stats["gene_3p"],
                "target_pair": [target_stats["gene_5p"], target_stats["gene_3p"]],
                "observed_count": target_stats["observed_count"],
                "expected_count": target_stats["expected_count"],
                "marginal_5p_count": target_stats["marginal_5p_count"],
                "marginal_3p_count": target_stats["marginal_3p_count"],
                "fold_enrichment": target_stats["fold_enrichment"],
                "odds_ratio": target_stats["odds_ratio"],
                "p_value": target_stats["p_value"],
                "is_significant": target_stats["is_significant"],
            })
        elif all_pair_stats:
            top = all_pair_stats[0]
            summary.update({
                "top_gene_5p": top["gene_5p"],
                "top_gene_3p": top["gene_3p"],
                "top_p_value": top["p_value"],
                "top_is_significant": top["is_significant"],
            })

        tables = {
            "pair_stats": all_pair_stats,
        }

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters=params,
            Summary=summary,
            Tables=tables,
            Warnings=warnings if warnings else None,
            Created_at=datetime.now(timezone.utc),
        )


JointPartnerMode = JointPartnerAlgorithm

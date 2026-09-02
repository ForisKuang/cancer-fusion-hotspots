"""Statistical tests for joint-partner-dependent fusion co-occurrence.

This module evaluates whether a specific gene pair (or pair orientation)
occurs significantly more often than would be predicted under an independence
null model based on marginal partner frequencies across the cohort.

# NOTE on Independence Null Simplification:
# The independence null model used here assumes that partner gene pairing occurs
# randomly in proportion to the marginal fusion frequency of each partner gene across
# the cohort: P(gene5, gene3) = P(gene5) * P(gene3).
#
# This is explicitly a first-pass simplification and stub for joint-partner dependency
# analysis:
# 1. Real biological joint-dependence involves tissue-specific chromatin architecture,
#    spatial genomic proximity (e.g. topological associating domains / TADs, chromatin loops),
#    and selective viability constraints that favor specific functional chimeras over
#    random pairings.
# 2. Breakpoint mechanisms (e.g. non-homologous end joining, fragile sites) create
#    non-uniform partner landscapes independent of oncogenic selection.
# 3. Sequencing panel and assay ascertainment biases often pre-select known partner
#    targets.
# Future iterations will extend this baseline with panel-aware and chromatin-aware nulls.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Optional

import numpy as np
from scipy.stats import binomtest, fisher_exact

from cfh.model.fusion_event import FusionEvent


def extract_event_genes(
    event: FusionEvent,
    directional: bool = True,
) -> Optional[tuple[str, str]]:
    """Extract partner genes from a FusionEvent.

    If directional is True:
        Attempts to extract (Five_prime_gene, Three_prime_gene).
        Falls back to (Site1_gene, Site2_gene).
    If directional is False:
        Sorts the extracted partner pair alphabetically.

    Returns None if fewer than two gene symbols can be identified.
    """
    g5: Optional[str] = None
    g3: Optional[str] = None

    if event.Five_prime_gene and event.Three_prime_gene:
        g5 = str(event.Five_prime_gene).strip().upper()
        g3 = str(event.Three_prime_gene).strip().upper()
    elif event.Site1_gene and event.Site2_gene:
        g5 = str(event.Site1_gene).strip().upper()
        g3 = str(event.Site2_gene).strip().upper()

    if not g5 or not g3:
        return None

    if not directional:
        return (min(g5, g3), max(g5, g3))
    return (g5, g3)


def compute_pair_enrichment(
    observed_count: int,
    total_events: int,
    marginal_5p_count: int,
    marginal_3p_count: int,
    method: str = "fisher",
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Test whether observed co-occurrence exceeds the independence null.

    Args:
        observed_count: Number of events with (gene5, gene3).
        total_events: Total informative fusion events in the cohort.
        marginal_5p_count: Total events where 5' partner is gene5.
        marginal_3p_count: Total events where 3' partner is gene3.
        method: Statistical test to perform ('fisher' for Fisher's exact test,
            'binomial' for binomial test against independence rate).
        alpha: Significance threshold (default: 0.05).

    Returns:
        dict containing observed, expected, fold enrichment, p-value, and test metadata.
    """
    if total_events <= 0:
        return {
            "observed_count": observed_count,
            "total_events": total_events,
            "marginal_5p_count": marginal_5p_count,
            "marginal_3p_count": marginal_3p_count,
            "expected_count": 0.0,
            "fold_enrichment": 0.0,
            "odds_ratio": None,
            "p_value": 1.0,
            "is_significant": False,
            "method": method,
            "alpha": alpha,
        }

    expected_count = (marginal_5p_count * marginal_3p_count) / total_events
    p_expected = (marginal_5p_count / total_events) * (marginal_3p_count / total_events)

    if expected_count > 0:
        fold_enrichment = float(observed_count / expected_count)
    elif observed_count > 0:
        fold_enrichment = float("inf")
    else:
        fold_enrichment = 0.0

    odds_ratio: Optional[float] = None
    chosen_method = method.lower().strip()

    if chosen_method == "fisher":
        a = observed_count
        b = max(0, marginal_5p_count - observed_count)
        c = max(0, marginal_3p_count - observed_count)
        d = max(0, total_events - marginal_5p_count - marginal_3p_count + observed_count)

        if a + b + c + d == 0:
            p_value = 1.0
        else:
            table = [[a, b], [c, d]]
            res = fisher_exact(table, alternative="greater")
            p_value = float(res.pvalue)
            stat = res.statistic
            if stat is not None and not np.isnan(stat) and not math.isnan(stat):
                odds_ratio = float(stat) if not math.isinf(stat) else float("inf")
    elif chosen_method == "binomial":
        if p_expected <= 0.0:
            p_value = 0.0 if observed_count > 0 else 1.0
        elif p_expected >= 1.0:
            p_value = 1.0
        else:
            res_b = binomtest(observed_count, total_events, p=p_expected, alternative="greater")
            p_value = float(res_b.pvalue)
    else:
        raise ValueError(f"Unsupported method {method!r}. Choose 'fisher' or 'binomial'.")

    return {
        "observed_count": observed_count,
        "total_events": total_events,
        "marginal_5p_count": marginal_5p_count,
        "marginal_3p_count": marginal_3p_count,
        "expected_count": float(expected_count),
        "fold_enrichment": fold_enrichment,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "is_significant": bool(p_value < alpha),
        "method": chosen_method,
        "alpha": alpha,
    }


def evaluate_gene_pair(
    events: list[FusionEvent],
    gene_5p: str,
    gene_3p: str,
    directional: bool = True,
    method: str = "fisher",
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Extract pairs from events and test a specific gene pair for co-occurrence enrichment."""
    target_5p = gene_5p.strip().upper()
    target_3p = gene_3p.strip().upper()
    if not directional:
        target_5p, target_3p = min(target_5p, target_3p), max(target_5p, target_3p)

    pairs: list[tuple[str, str]] = []
    for ev in events:
        p = extract_event_genes(ev, directional=directional)
        if p is not None:
            pairs.append(p)

    total_events = len(pairs)
    observed_count = sum(1 for p in pairs if p == (target_5p, target_3p))
    marginal_5p_count = sum(1 for p in pairs if p[0] == target_5p)
    marginal_3p_count = sum(1 for p in pairs if p[1] == target_3p)

    stats = compute_pair_enrichment(
        observed_count=observed_count,
        total_events=total_events,
        marginal_5p_count=marginal_5p_count,
        marginal_3p_count=marginal_3p_count,
        method=method,
        alpha=alpha,
    )
    stats["gene_5p"] = target_5p
    stats["gene_3p"] = target_3p
    return stats


def evaluate_all_pairs(
    events: list[FusionEvent],
    directional: bool = True,
    method: str = "fisher",
    alpha: float = 0.05,
    min_count: int = 1,
) -> list[dict[str, Any]]:
    """Extract all unique gene pairs from events and compute co-occurrence enrichment for each."""
    pairs: list[tuple[str, str]] = []
    for ev in events:
        p = extract_event_genes(ev, directional=directional)
        if p is not None:
            pairs.append(p)

    total_events = len(pairs)
    if total_events == 0:
        return []

    pair_counts = Counter(pairs)
    marginal_5p = Counter(p[0] for p in pairs)
    marginal_3p = Counter(p[1] for p in pairs)

    results: list[dict[str, Any]] = []
    for (g5, g3), obs in pair_counts.items():
        if obs < min_count:
            continue
        st = compute_pair_enrichment(
            observed_count=obs,
            total_events=total_events,
            marginal_5p_count=marginal_5p[g5],
            marginal_3p_count=marginal_3p[g3],
            method=method,
            alpha=alpha,
        )
        st["gene_5p"] = g5
        st["gene_3p"] = g3
        results.append(st)

    # Sort by p_value ascending, then observed_count descending
    results.sort(key=lambda r: (r["p_value"], -r["observed_count"]))
    return results

"""Gene-agnostic, frequency/recurrence-only cutpoint (boundary) detection.

No external oncogenicity label and no OncoKB dependency: the outcome used
here is domain-retention status (retained vs. lost/disrupted), already
computed by the domain-retention algorithm and stored per-event via
``FusionFeature.Domain_retention_flags``. The scanning axis is each event's
already-mapped breakpoint protein position
(``FusionFeature.Junction_position_aa``); recurrence is simply how many
independent events land at or near a given position. No new network call is
introduced: an optional ``GenomeNexusClient`` may be passed in ``params``
purely to compare the inferred cutpoint against Pfam domain boundaries the
pipeline already fetches elsewhere, and the algorithm runs fully offline
without one.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from cfh.algorithms.base import Algorithm
from cfh.algorithms.registry import register
from cfh.genes.registry import GeneConfig
from cfh.mapping.genome_nexus_source import GenomeNexusClient, parse_canonical_transcript
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature
from cfh.stats.breakpoint_tests import gene_breakpoint_domain_status_records
from cfh.stats.cutpoint_scan import detect_cutpoint

ALGORITHM_NAME = "cutpoint_detection"
ALGORITHM_VERSION = "0.1.0"


def _nearest_boundary_comparison(cutpoint: int, boundaries: list[int]) -> Optional[dict]:
    if not boundaries:
        return None
    nearest = min(boundaries, key=lambda boundary: abs(boundary - cutpoint))
    return {
        "nearest_known_domain_boundary_aa": nearest,
        "distance_aa": abs(nearest - cutpoint),
    }


def _known_domain_boundaries(gene_config: GeneConfig, params: dict) -> list[int]:
    """Known Pfam/domain boundary positions to validate the inferred cutpoint against.

    Prefers an explicit ``params["domain_boundaries"]`` list (e.g. sourced
    from a committed fixture) so the benchmark stays reproducible offline.
    Falls back to an optional ``params["genome_nexus_client"]`` -- the same
    client type ``domain_retention`` already accepts -- to fetch Pfam
    domain start/end positions for the gene's canonical transcript. Neither
    is required; with neither supplied the comparison is simply omitted.
    """
    explicit = params.get("domain_boundaries")
    if explicit:
        return sorted({int(boundary) for boundary in explicit})

    client: Optional[GenomeNexusClient] = params.get("genome_nexus_client")
    if client is None or not gene_config.gene_symbol:
        return []

    payload = client.fetch_canonical_transcript(gene_config.gene_symbol)
    transcript = parse_canonical_transcript(payload)
    boundaries: set[int] = set()
    for domain in transcript.pfam_domains:
        boundaries.add(domain.start_aa)
        boundaries.add(domain.end_aa)
    return sorted(boundaries)


@register(ALGORITHM_NAME)
class CutpointDetectionAlgorithm(Algorithm):
    """Scan a gene's observed breakpoints for the protein-position cutpoint
    that best separates domain-retained from lost/disrupted events, with a
    permutation-corrected empirical p-value for having scanned many
    candidate cutpoints.

    Expected ``params`` keys (all optional):
        seed (int): permutation RNG seed, default 42.
        n_permutations (int): number of label permutations, default 10000.
        domain_boundaries (list[int]): known domain boundary aa positions to
            compare the inferred cutpoint against.
        genome_nexus_client (GenomeNexusClient): used only as a fallback
            source for ``domain_boundaries`` when that param is omitted.
    """

    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        params = params or {}
        seed = params.get("seed", 42)
        n_permutations = params.get("n_permutations", 10_000)

        records = gene_breakpoint_domain_status_records(events, features, gene_config)
        positions = [position for position, _ in records]
        statuses = [status for _, status in records]

        scan_result = detect_cutpoint(
            positions,
            statuses,
            seed=seed,
            n_permutations=n_permutations,
        )

        boundary_comparison = None
        warnings: list[str] = []
        if scan_result["determinable"]:
            boundaries = _known_domain_boundaries(gene_config, params)
            boundary_comparison = _nearest_boundary_comparison(
                scan_result["best_cutpoint"], boundaries
            )
        else:
            warnings.append(scan_result["reason"])

        return AlgorithmResult(
            Algorithm=ALGORITHM_NAME,
            Algorithm_version=ALGORITHM_VERSION,
            Parameters={"seed": seed, "n_permutations": n_permutations},
            Summary={
                "determinable": scan_result["determinable"],
                "reason": scan_result["reason"],
                "n_events_analyzed": len(records),
                "inferred_cutpoint_aa": scan_result["best_cutpoint"],
                "observed_statistic_neg_log10_p": scan_result["observed_statistic"],
                "observed_p_value": scan_result["observed_p_value"],
                "observed_odds_ratio": scan_result["observed_odds_ratio"],
                "corrected_p_value": scan_result["corrected_p_value"],
                "known_domain_boundary_comparison": boundary_comparison,
            },
            Tables={"cutpoint_scan": scan_result["scan"]},
            Warnings=warnings,
            Created_at=datetime.now(timezone.utc),
        )

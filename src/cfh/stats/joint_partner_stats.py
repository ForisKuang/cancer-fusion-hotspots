"""Pair co-occurrence enrichment for joint-partner fusion analysis."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from scipy.stats import fisher_exact

from cfh.model.fusion_event import FusionEvent


@dataclass(frozen=True)
class PairEnrichment:
    """Observed and null-model statistics for one ordered fusion-gene pair."""

    gene5: str
    gene3: str
    eligible_event_count: int
    observed_count: int
    expected_count: float
    p_value: float
    odds_ratio: float

    def as_dict(self) -> dict[str, float | int | str]:
        """Return a JSON-serializable representation for ``AlgorithmResult.Tables``."""
        return asdict(self)


def _event_pair(event: FusionEvent) -> tuple[str, str] | None:
    """Return the oriented pair, falling back to the source breakpoint-site order."""
    if event.Five_prime_gene is not None and event.Three_prime_gene is not None:
        return event.Five_prime_gene, event.Three_prime_gene
    if event.Site1_gene is not None and event.Site2_gene is not None:
        return event.Site1_gene, event.Site2_gene
    return None


def calculate_pair_enrichment(events: list[FusionEvent], gene5: str, gene3: str) -> PairEnrichment:
    """Test whether an ordered pair exceeds its marginal-independence expectation.

    Events without a complete gene pair are excluded. The Fisher exact test is
    one-sided (``greater``), since this mode looks only for overrepresentation.
    """
    pairs = [pair for event in events if (pair := _event_pair(event)) is not None]
    total = len(pairs)
    pair_count = sum(pair == (gene5, gene3) for pair in pairs)
    gene5_count = sum(pair[0] == gene5 for pair in pairs)
    gene3_count = sum(pair[1] == gene3 for pair in pairs)

    expected_count = (gene5_count * gene3_count / total) if total else 0.0
    if total == 0:
        return PairEnrichment(
            gene5=gene5,
            gene3=gene3,
            eligible_event_count=0,
            observed_count=0,
            expected_count=expected_count,
            p_value=1.0,
            odds_ratio=0.0,
        )

    # This marginal-independence null is intentionally a simplification:
    # real biological joint-dependence is more complex than independent
    # partner selection in a cohort.
    contingency_table = [
        [pair_count, gene5_count - pair_count],
        [gene3_count - pair_count, total - gene5_count - gene3_count + pair_count],
    ]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")
    return PairEnrichment(
        gene5=gene5,
        gene3=gene3,
        eligible_event_count=total,
        observed_count=pair_count,
        expected_count=expected_count,
        p_value=float(p_value),
        odds_ratio=float(odds_ratio),
    )

"""Map a FusionEvent + GeneConfig onto a FusionFeature, including per-domain
retention status.

Generic by construction: every gene-specific fact (which domains exist,
what they're called, which one is which) comes from the ``GeneConfig`` and
the fetched domain source, never from a literal in this module.
"""

from __future__ import annotations

from cfh.genes.registry import GeneConfig, KeyDomain
from cfh.mapping.domain_source import ProteinDomain, UniProtDomainSource
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature

_default_domain_source: UniProtDomainSource | None = None


def get_default_domain_source() -> UniProtDomainSource:
    """Return the process-wide default :class:`UniProtDomainSource`.

    ``map_event`` uses this whenever the caller doesn't pass an explicit
    ``domain_source``, so a batch of events processed through the real
    call path (no manually-shared instance required) still shares one
    cache and therefore one HTTP call per accession, not one per event.
    """
    global _default_domain_source
    if _default_domain_source is None:
        _default_domain_source = UniProtDomainSource()
    return _default_domain_source


def reset_default_domain_source() -> None:
    """Drop the cached default domain source (mainly for test isolation)."""
    global _default_domain_source
    _default_domain_source = None


def _normalize_domain_name(name: str) -> str:
    normalized = name.lower().strip()
    for suffix in (" domain", " region"):
        if normalized.endswith(suffix):
            normalized = normalized[: -len(suffix)]
    return normalized.strip()


def _find_matching_domain(
    domains: list[ProteinDomain], config_domain: KeyDomain
) -> ProteinDomain | None:
    target = _normalize_domain_name(config_domain.name)
    for domain in domains:
        if config_domain.accession and config_domain.accession in {
            domain.accession,
            domain.name,
        }:
            return domain
        candidate = _normalize_domain_name(domain.name)
        if target == candidate or target in candidate or candidate in target:
            return domain
    return None


def classify_domain_retention(
    domain_start: int | None,
    domain_end: int | None,
    breakpoint_aa: int | None,
    role: str | None,
) -> str:
    """Classify a domain as retained/disrupted/lost given a breakpoint and role.

    ``role`` is ``"five_prime"`` (the retained fragment runs from the
    protein start up to the breakpoint) or ``"three_prime"`` (the retained
    fragment runs from the breakpoint to the protein end).
    """
    if breakpoint_aa is None or domain_start is None or domain_end is None:
        return "unknown"

    if role == "five_prime":
        retained_start, retained_end = 0, breakpoint_aa
    elif role == "three_prime":
        retained_start, retained_end = breakpoint_aa, float("inf")
    else:
        return "unknown"

    if domain_start >= retained_start and domain_end <= retained_end:
        return "retained"
    if domain_start > retained_end or domain_end < retained_start:
        return "lost"
    return "disrupted"


def map_event(
    event: FusionEvent,
    gene_config: GeneConfig,
    *,
    role: str,
    junction_position_aa: int | None,
    domain_source: UniProtDomainSource | None = None,
) -> FusionFeature:
    """Build a FusionFeature for ``gene_config``'s gene in ``event``.

    ``role`` and ``junction_position_aa`` describe this gene's side of the
    breakpoint in protein-amino-acid coordinates (from transcript/exon
    mapping, out of scope here).
    """
    domain_source = domain_source or get_default_domain_source()
    domains = domain_source.fetch(gene_config.protein_id)

    retention_flags: dict[str, str] = {}
    retained_domains: list[str] = []
    lost_domains: list[str] = []
    disrupted_domains: list[str] = []

    for key_domain in [*gene_config.key_domains, *gene_config.disruption_required_domains]:
        matched = _find_matching_domain(domains, key_domain)
        flag_key = key_domain.key or _normalize_domain_name(key_domain.name)
        status = classify_domain_retention(
            matched.start_aa if matched else None,
            matched.end_aa if matched else None,
            junction_position_aa,
            role,
        )
        retention_flags[flag_key] = status
        if status == "retained":
            retained_domains.append(key_domain.name)
        elif status == "lost":
            lost_domains.append(key_domain.name)
        elif status == "disrupted":
            disrupted_domains.append(key_domain.name)

    return FusionFeature(
        Event_id=event.Event_id,
        Gene=gene_config.gene_symbol,
        Role=role,
        Transcript_id=gene_config.canonical_transcript_id,
        Junction_position_aa=junction_position_aa,
        Retained_domains=retained_domains,
        Lost_domains=lost_domains,
        Disrupted_domains=disrupted_domains,
        Domain_retention_flags=retention_flags,
    )

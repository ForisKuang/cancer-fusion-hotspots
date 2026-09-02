"""Resolve a fusion breakpoint to a transcript/exon position.

Primary source: exon numbers already present in the cBioPortal SV
annotation text. Fallback: Ensembl's REST overlap API for the protein's
features, used when the annotation doesn't mention an exon at all.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

import requests

from cfh.genes.registry import GeneConfig

_EXON_PATTERN = re.compile(r"exon\s*(\d+)", re.IGNORECASE)


class EnsemblFallbackUnavailable(RuntimeError):
    """Raised when the Ensembl-fallback mapping cannot even be attempted.

    This is the case when the SV annotation has no exon number AND the
    caller has no ``ensembl_protein_id`` to query Ensembl with. Callers
    should treat this as "fallback unavailable" and must not confuse it
    with "queried Ensembl and it had nothing to say" (which is a real,
    non-exceptional result -- see :func:`resolve_transcript_mapping`).
    """


@dataclass
class TranscriptMapping:
    transcript_id: str
    breakpoint_exon: int | None
    source: str  # "annotation" | "ensembl_fallback"
    ensembl_features: list[dict] | None = None
    breakpoint_protein_features: list[dict] | None = None
    """Ensembl protein features (e.g. domains) whose range contains the
    supplied breakpoint amino-acid coordinate, when one was given. This is
    the real, computed result of the Ensembl fallback: the
    ``overlap/translation`` endpoint exposes protein-feature boundaries,
    not an exon-to-protein-coordinate map, so ``breakpoint_exon`` genuinely
    stays ``None`` on this path -- that is a data-source limitation, not a
    silently-dropped answer.
    """


class EnsemblClient:
    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, session: "requests.Session | None" = None):
        self.session = session or requests.Session()

    def fetch_protein_features(self, ensembl_protein_id: str) -> list[dict]:
        url = f"{self.BASE_URL}/overlap/translation/{ensembl_protein_id}"
        response = self.session.get(
            url,
            params={"feature": "protein_feature"},
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        response.raise_for_status()
        return response.json()


def extract_exon_from_annotation(annotation: str | None) -> int | None:
    """Pull an explicit exon number out of free-text SV annotation, if present."""
    if not annotation:
        return None
    match = _EXON_PATTERN.search(annotation)
    return int(match.group(1)) if match else None


def _features_overlapping(features: list[dict], breakpoint_aa: int) -> list[dict]:
    overlapping = []
    for feature in features:
        start, end = feature.get("start"), feature.get("end")
        if start is None or end is None:
            continue
        if start <= breakpoint_aa <= end:
            overlapping.append(feature)
    return overlapping


def resolve_transcript_mapping(
    annotation: str | None,
    gene_config: GeneConfig,
    *,
    breakpoint_aa: int | None = None,
    ensembl_protein_id: str | None = None,
    ensembl_client: "EnsemblClient | Any | None" = None,
) -> TranscriptMapping:
    """Resolve breakpoint exon info, preferring the SV annotation text.

    When the annotation lacks an explicit exon number, this falls back to
    fetching Ensembl protein features for ``ensembl_protein_id`` and, when a
    ``breakpoint_aa`` protein-coordinate is supplied, actually uses those
    fetched features to report which ones overlap the breakpoint
    (``breakpoint_protein_features``) -- a real result computed from the
    fallback data source, not a silent ``None``.

    Raises:
        EnsemblFallbackUnavailable: the annotation has no exon number and
            no ``ensembl_protein_id`` was supplied, so the fallback cannot
            be attempted at all.
    """
    exon = extract_exon_from_annotation(annotation)
    if exon is not None:
        return TranscriptMapping(
            transcript_id=gene_config.canonical_transcript_id,
            breakpoint_exon=exon,
            source="annotation",
        )

    if not ensembl_protein_id:
        raise EnsemblFallbackUnavailable(
            "annotation has no exon number and no ensembl_protein_id was "
            "supplied, so the Ensembl fallback cannot be attempted"
        )

    client = ensembl_client or EnsemblClient()
    features = client.fetch_protein_features(ensembl_protein_id)

    breakpoint_protein_features = (
        _features_overlapping(features, breakpoint_aa) if breakpoint_aa is not None else None
    )

    return TranscriptMapping(
        transcript_id=gene_config.canonical_transcript_id,
        breakpoint_exon=None,
        source="ensembl_fallback",
        ensembl_features=features,
        breakpoint_protein_features=breakpoint_protein_features,
    )

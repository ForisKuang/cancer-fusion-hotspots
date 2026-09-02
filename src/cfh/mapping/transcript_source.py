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


@dataclass
class TranscriptMapping:
    transcript_id: str
    breakpoint_exon: int | None
    source: str  # "annotation" | "ensembl_fallback"
    ensembl_features: list[dict] | None = None


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


def resolve_transcript_mapping(
    annotation: str | None,
    gene_config: GeneConfig,
    *,
    ensembl_protein_id: str | None = None,
    ensembl_client: "EnsemblClient | Any | None" = None,
) -> TranscriptMapping:
    """Resolve breakpoint exon info, preferring the SV annotation text."""
    exon = extract_exon_from_annotation(annotation)
    if exon is not None:
        return TranscriptMapping(
            transcript_id=gene_config.canonical_transcript_id,
            breakpoint_exon=exon,
            source="annotation",
        )

    features: list[dict] = []
    if ensembl_protein_id:
        client = ensembl_client or EnsemblClient()
        features = client.fetch_protein_features(ensembl_protein_id)

    return TranscriptMapping(
        transcript_id=gene_config.canonical_transcript_id,
        breakpoint_exon=None,
        source="ensembl_fallback",
        ensembl_features=features,
    )

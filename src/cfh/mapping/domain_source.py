"""Protein domain sources: UniProt (primary) and InterPro (cross-check/fallback).

Both accept any UniProt/InterPro accession generically; no gene-specific
logic lives here. An in-memory (and optional on-disk) cache means a given
accession triggers at most one real HTTP request per process/cache
directory, which keeps the pipeline reproducible offline once populated.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import requests


@dataclass
class ProteinDomain:
    name: str
    start_aa: int
    end_aa: int
    source: str  # "uniprot" | "interpro" | "genome_nexus"
    accession: str | None = None


_DOMAIN_FEATURE_TYPES = {"domain", "region"}


class UniProtDomainSource:
    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(
        self,
        session: "requests.Session | None" = None,
        cache_dir: str | Path | None = None,
    ):
        self.session = session or requests.Session()
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self._cache: dict[str, list[ProteinDomain]] = {}

    def parse(self, payload: dict) -> list[ProteinDomain]:
        """Parse a UniProt ``/uniprotkb/{accession}.json`` payload into domains."""
        domains: list[ProteinDomain] = []
        for feature in payload.get("features", []):
            if feature.get("type", "").lower() not in _DOMAIN_FEATURE_TYPES:
                continue
            location = feature.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            description = feature.get("description")
            if description is None or start is None or end is None:
                continue
            domains.append(
                ProteinDomain(
                    name=description,
                    start_aa=int(start),
                    end_aa=int(end),
                    source="uniprot",
                )
            )
        return domains

    def _cache_file(self, accession: str) -> Path | None:
        return self.cache_dir / f"{accession}.json" if self.cache_dir else None

    def fetch(self, accession: str) -> list[ProteinDomain]:
        """Fetch (and cache) domains for a UniProt accession, e.g. ``"P15056"``."""
        if accession in self._cache:
            return self._cache[accession]

        cache_file = self._cache_file(accession)
        if cache_file and cache_file.exists():
            payload = json.loads(cache_file.read_text())
        else:
            url = f"{self.BASE_URL}/{accession}.json"
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            payload = response.json()
            if cache_file:
                cache_file.parent.mkdir(parents=True, exist_ok=True)
                cache_file.write_text(json.dumps(payload))

        domains = self.parse(payload)
        self._cache[accession] = domains
        return domains


class InterProDomainSource:
    """Cross-check / fallback domain source using InterPro's Pfam-in-UniProt API."""

    BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot"

    def __init__(self, session: "requests.Session | None" = None):
        self.session = session or requests.Session()
        self._cache: dict[str, list[ProteinDomain]] = {}

    def parse(self, payload: dict) -> list[ProteinDomain]:
        domains: list[ProteinDomain] = []
        for result in payload.get("results", []):
            metadata = result.get("metadata", {})
            name = metadata.get("name")
            for entry_protein in result.get("proteins", []):
                for location in entry_protein.get("entry_protein_locations", []):
                    for fragment in location.get("fragments", []):
                        start, end = fragment.get("start"), fragment.get("end")
                        if name is None or start is None or end is None:
                            continue
                        domains.append(
                            ProteinDomain(
                                name=name,
                                start_aa=int(start),
                                end_aa=int(end),
                                source="interpro",
                            )
                        )
        return domains

    def fetch(self, accession: str) -> list[ProteinDomain]:
        if accession in self._cache:
            return self._cache[accession]
        url = f"{self.BASE_URL}/{accession}/"
        response = self.session.get(url, timeout=30)
        response.raise_for_status()
        domains = self.parse(response.json())
        self._cache[accession] = domains
        return domains

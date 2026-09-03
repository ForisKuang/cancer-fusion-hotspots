"""Per-gene biological configuration registry.

Gene-specific facts (canonical transcript, protein accession, domain
boundaries, ...) live in YAML files under ``genes/configs/`` and are loaded
into a validated :class:`GeneConfig`. Generic pipeline code must never
hardcode a gene's biology directly; it should always receive a
``GeneConfig`` instance instead.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, ConfigDict, model_validator

CONFIGS_DIR = Path(__file__).parent / "configs"


class KeyDomain(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    source: str
    key: Optional[str] = None
    accession: Optional[str] = None
    """Source-native identifier, such as a Pfam accession from Genome Nexus."""


class GeneConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene_symbol: Optional[str] = None
    gene_pair: Optional[tuple[str, str]] = None
    canonical_transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    key_domains: list[KeyDomain] = []
    autoinhibitory_domains: list[str] = []
    expected_retained_exon_hint: Optional[str] = None
    analysis_modes: list[str] = []
    entrez_gene_id: Optional[int] = None
    """NCBI Entrez gene id, e.g. for cBioPortal structural-variant API queries."""

    @model_validator(mode="after")
    def _validate_config_target(self) -> "GeneConfig":
        """Require either a single gene or an ordered fusion-gene pair."""
        if self.gene_symbol is None and self.gene_pair is None:
            raise ValueError("GeneConfig requires gene_symbol or gene_pair")
        if self.gene_pair is None and (
            self.canonical_transcript_id is None or self.protein_id is None
        ):
            raise ValueError(
                "single-gene GeneConfig requires canonical_transcript_id and protein_id"
            )
        return self


def _config_path(gene_symbol: str) -> Path:
    return CONFIGS_DIR / f"{gene_symbol.lower()}.yaml"


def load_gene_config(gene_symbol: str) -> GeneConfig:
    """Load and validate a gene's YAML config by symbol (e.g. ``"braf"``)."""
    path = _config_path(gene_symbol)
    if not path.exists():
        raise FileNotFoundError(f"No gene config found for {gene_symbol!r} at {path}")
    with path.open() as fh:
        raw = yaml.safe_load(fh)
    return GeneConfig.model_validate(raw)


def available_genes() -> list[str]:
    """List gene symbols with a registered config."""
    return sorted(p.stem.upper() for p in CONFIGS_DIR.glob("*.yaml"))

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
from pydantic import BaseModel, ConfigDict

CONFIGS_DIR = Path(__file__).parent / "configs"


class KeyDomain(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    source: str
    key: Optional[str] = None


class GeneConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene_symbol: str
    canonical_transcript_id: Optional[str] = None
    protein_id: Optional[str] = None
    key_domains: list[KeyDomain] = []
    autoinhibitory_domains: list[str] = []
    expected_retained_exon_hint: Optional[str] = None
    analysis_modes: list[str] = []
    entrez_gene_id: Optional[int] = None
    """NCBI Entrez gene id, e.g. for cBioPortal structural-variant API queries."""
    gene_pair: Optional[list[str]] = None
    partner_5p: Optional[str] = None
    partner_3p: Optional[str] = None
    description: Optional[str] = None


def _config_path(gene_symbol: str) -> Path:
    target = CONFIGS_DIR / f"{gene_symbol.lower()}.yaml"
    if not target.exists():
        hyphenated = CONFIGS_DIR / f"{gene_symbol.lower().replace('_', '-')}.yaml"
        if hyphenated.exists():
            return hyphenated
    return target


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

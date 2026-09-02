"""Genome Nexus integration: canonical-transcript/Pfam-domain lookup by gene
symbol, and genomic-breakpoint-to-protein-position mapping using a
transcript's exon structure.

Generic by construction: every function here takes a gene symbol,
accession, or exon list supplied by the caller (typically from a
``GeneConfig`` or a fetched Genome Nexus payload) -- nothing gene-specific
is hardcoded.

Genome Nexus's own variant-annotation endpoints are for point mutations
expressed as HGVS/genomic-change strings, not structural-variant fusion
breakpoints -- there is no endpoint that maps a fusion breakpoint to a
protein position directly. This module does that arithmetic itself (see
:func:`map_genomic_breakpoint_to_protein_position`) from the exon
coordinates Genome Nexus does provide.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

import requests

from cfh.mapping.domain_source import ProteinDomain

DEFAULT_BASE_URL = "https://www.genomenexus.org"
DEFAULT_ISOFORM_OVERRIDE_SOURCE = "mskcc"

_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


class GenomeNexusGeneNotFound(RuntimeError):
    """Raised when Genome Nexus has no canonical-transcript mapping at all
    for a requested gene symbol (a 404 from the canonical-transcript
    endpoint). Callers should treat this as "try the next domain/transcript
    source" (e.g. UniProt-direct), not as a transient failure.
    """


@dataclass
class PfamDomain:
    pfam_id: str
    start_aa: int
    end_aa: int


@dataclass
class ExonRecord:
    exon_id: str
    start: int
    end: int
    rank: int
    strand: int


@dataclass
class CanonicalTranscript:
    transcript_id: str
    refseq_mrna_id: str | None
    protein_id: str | None
    protein_length: int | None
    uniprot_id: str | None
    pfam_domains: list[PfamDomain]
    exons: list[ExonRecord]


@dataclass
class BreakpointProteinPosition:
    protein_position: int
    exon_rank: int
    is_intronic: bool
    """True if the breakpoint fell in an intron and had to be clamped to
    the nearest exon boundary rather than landing inside a coding exon.
    """


class GenomeNexusClient:
    """Thin REST client with basic retry/backoff and a per-process cache.

    The public Genome Nexus instance documents no rate limit and requires
    no auth, but repeatedly querying the same gene/transcript in a batch
    pipeline is still wasteful and impolite, so results are cached by
    request key for the life of the client instance.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        session: "requests.Session | None" = None,
        max_retries: int = 3,
        backoff_seconds: float = 0.5,
    ):
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self._canonical_transcript_cache: dict[str, dict] = {}
        self._transcript_cache: dict[str, dict] = {}

    def _get(self, path: str, *, params: dict | None = None) -> requests.Response:
        url = f"{self.base_url}{path}"
        attempt = 0
        while True:
            response = self.session.get(url, params=params, timeout=30)
            if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= self.max_retries:
                return response
            time.sleep(self.backoff_seconds * (2**attempt))
            attempt += 1

    def fetch_canonical_transcript(
        self, hugo_symbol: str, *, isoform_override_source: str = DEFAULT_ISOFORM_OVERRIDE_SOURCE
    ) -> dict:
        """GET ``/ensembl/canonical-transcript/hgnc/{hugoSymbol}``, cached by symbol."""
        cache_key = f"{hugo_symbol}:{isoform_override_source}"
        if cache_key in self._canonical_transcript_cache:
            return self._canonical_transcript_cache[cache_key]

        response = self._get(
            f"/ensembl/canonical-transcript/hgnc/{hugo_symbol}",
            params={"isoformOverrideSource": isoform_override_source},
        )
        if response.status_code == 404:
            raise GenomeNexusGeneNotFound(
                f"Genome Nexus has no canonical-transcript mapping for gene symbol {hugo_symbol!r}"
            )
        response.raise_for_status()
        payload = response.json()
        self._canonical_transcript_cache[cache_key] = payload
        return payload

    def fetch_transcript(self, transcript_id: str) -> dict:
        """GET ``/ensembl/transcript/{transcriptId}``, cached by transcript id."""
        if transcript_id in self._transcript_cache:
            return self._transcript_cache[transcript_id]

        response = self._get(f"/ensembl/transcript/{transcript_id}")
        response.raise_for_status()
        payload = response.json()
        self._transcript_cache[transcript_id] = payload
        return payload


def parse_canonical_transcript(payload: dict) -> CanonicalTranscript:
    """Parse a canonical-transcript (or transcript-by-id) Genome Nexus payload."""
    pfam_domains = [
        PfamDomain(
            pfam_id=domain["pfamDomainId"],
            start_aa=domain["pfamDomainStart"],
            end_aa=domain["pfamDomainEnd"],
        )
        for domain in payload.get("pfamDomains", [])
    ]
    exons = [
        ExonRecord(
            exon_id=exon["exonId"],
            start=exon["exonStart"],
            end=exon["exonEnd"],
            rank=exon["rank"],
            strand=exon["strand"],
        )
        for exon in payload.get("exons", [])
    ]
    return CanonicalTranscript(
        transcript_id=payload["transcriptId"],
        refseq_mrna_id=payload.get("refseqMrnaId"),
        protein_id=payload.get("proteinId"),
        protein_length=payload.get("proteinLength"),
        uniprot_id=payload.get("uniprotId"),
        pfam_domains=pfam_domains,
        exons=exons,
    )


def pfam_domains_to_protein_domains(canonical: CanonicalTranscript) -> list[ProteinDomain]:
    """Adapt Genome Nexus Pfam domains to the shared :class:`ProteinDomain` shape."""
    return [
        ProteinDomain(
            name=domain.pfam_id,
            start_aa=domain.start_aa,
            end_aa=domain.end_aa,
            source="genome_nexus",
        )
        for domain in canonical.pfam_domains
    ]


def resolve_domains(
    gene_symbol: str,
    uniprot_accession: str | None,
    *,
    genome_nexus_client: "GenomeNexusClient | None" = None,
    uniprot_source=None,
) -> list[ProteinDomain]:
    """Domain-annotation source chain: Genome Nexus first, UniProt-direct fallback.

    Genome Nexus's canonical-transcript endpoint is queried by gene symbol
    and needs no accession lookup step, so it's tried first. If it has no
    mapping for this gene at all (:class:`GenomeNexusGeneNotFound`), this
    falls back to a direct UniProt lookup by accession rather than raising,
    so a gene missing from Genome Nexus still gets domain annotations.
    """
    client = genome_nexus_client or GenomeNexusClient()
    try:
        payload = client.fetch_canonical_transcript(gene_symbol)
    except GenomeNexusGeneNotFound:
        payload = None

    if payload is not None:
        canonical = parse_canonical_transcript(payload)
        domains = pfam_domains_to_protein_domains(canonical)
        if domains:
            return domains

    if uniprot_accession is None:
        return []

    from cfh.mapping.domain_source import UniProtDomainSource

    source = uniprot_source or UniProtDomainSource()
    return source.fetch(uniprot_accession)


def _coding_interval(
    exon: ExonRecord, cds_min_genomic: int | None, cds_max_genomic: int | None
) -> tuple[int, int] | None:
    start, end = exon.start, exon.end
    if cds_min_genomic is not None:
        start = max(start, cds_min_genomic)
    if cds_max_genomic is not None:
        end = min(end, cds_max_genomic)
    if start > end:
        return None  # exon is entirely outside the CDS (e.g. a pure-UTR exon)
    return start, end


def map_genomic_breakpoint_to_protein_position(
    exons: list[ExonRecord],
    breakpoint_genomic: int,
    strand: int,
    *,
    cds_min_genomic: int | None = None,
    cds_max_genomic: int | None = None,
) -> BreakpointProteinPosition:
    """Estimate the protein residue a genomic breakpoint falls in.

    This is standard breakpoint-to-CDS-position arithmetic, not something
    Genome Nexus (or any variant-annotation API) computes for a fusion
    breakpoint directly:

    1. Order exons 5'->3' along the mature transcript by ``rank`` ascending
       (transcript order, independent of genomic strand).
    2. Optionally clip each exon to ``[cds_min_genomic, cds_max_genomic]``
       so a first/last exon that partly overlaps the 5'/3' UTR only
       contributes its coding portion (an exon entirely outside this
       range contributes nothing).
    3. Find the exon whose clipped genomic range contains the breakpoint.
       If none does (an intronic breakpoint), use the nearest exon by
       genomic distance and clamp the breakpoint to that exon's nearest
       boundary -- an approximation, flagged via ``is_intronic``.
    4. Convert the breakpoint's position within that exon to a
       CDS-relative nucleotide offset: on the plus strand this counts up
       from the exon's clipped start; on the minus strand -- transcribed
       high-to-low in genomic coordinate space -- it counts down from the
       exon's clipped end.
    5. The 1-indexed CDS nucleotide position is the sum of the coding
       lengths of every exon ranked before this one, plus that offset,
       plus one.
    6. The protein residue number is that nucleotide position divided by
       3, rounded up (any nucleotide within a codon's 3 bases maps to
       that codon's amino acid).
    """
    ordered = sorted(exons, key=lambda e: e.rank)
    coding_intervals = {
        exon.rank: _coding_interval(exon, cds_min_genomic, cds_max_genomic) for exon in ordered
    }

    containing = None
    for exon in ordered:
        interval = coding_intervals[exon.rank]
        if interval is None:
            continue
        start, end = interval
        if start <= breakpoint_genomic <= end:
            containing = exon
            break

    is_intronic = containing is None
    if containing is None:

        def _distance(exon: ExonRecord) -> float:
            interval = coding_intervals[exon.rank]
            if interval is None:
                return math.inf
            start, end = interval
            if breakpoint_genomic < start:
                return start - breakpoint_genomic
            if breakpoint_genomic > end:
                return breakpoint_genomic - end
            return 0

        containing = min(ordered, key=_distance)
        start, end = coding_intervals[containing.rank]
        clamped_breakpoint = min(max(breakpoint_genomic, start), end)
    else:
        start, end = coding_intervals[containing.rank]
        clamped_breakpoint = breakpoint_genomic

    offset = (clamped_breakpoint - start) if strand == 1 else (end - clamped_breakpoint)

    preceding_length = sum(
        (coding_intervals[exon.rank][1] - coding_intervals[exon.rank][0] + 1)
        for exon in ordered
        if exon.rank < containing.rank and coding_intervals[exon.rank] is not None
    )

    cds_nt_position = preceding_length + offset + 1
    protein_position = math.ceil(cds_nt_position / 3)

    return BreakpointProteinPosition(
        protein_position=protein_position,
        exon_rank=containing.rank,
        is_intronic=is_intronic,
    )

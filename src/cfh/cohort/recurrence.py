"""Cohort-wide structural-variant gene recurrence and recurrence gating.

A genome-wide scan starts from cBioPortal's ``/structuralvariant-genes/fetch``
endpoint, which returns every gene that has at least one structural-variant
record anywhere in a study -- for ``msk_impact_50k_2026`` that is thousands
of genes, most seen in only a single patient (noise, not a recurrent
hotspot candidate). Recurrence gating filters that list down to genes
recurrent enough to be worth running the full per-gene algorithm suite on,
while always reporting the pre-filter total so ungated genes are never
silently dropped from view.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from cfh.ingestion import cbioportal_api

DEFAULT_MIN_DISTINCT_PATIENTS = 5


@dataclass(frozen=True)
class GeneRecurrence:
    """One gene's cohort-wide structural-variant recurrence."""

    hugo_gene_symbol: str
    entrez_gene_id: int | None
    distinct_patient_count: int
    """``numberOfAlteredCases`` from cBioPortal: the number of distinct
    patients with at least one structural-variant record touching this
    gene -- the recurrence signal the gate filters on."""
    total_sv_count: int
    """``totalCount``: the raw number of SV records for this gene, which can
    exceed ``distinct_patient_count`` (a patient with multiple SV records
    for the same gene) and is not itself the gating criterion."""


def fetch_cohort_gene_recurrence(
    study_id: str,
    *,
    base_url: str = cbioportal_api.DEFAULT_BASE_URL,
    session=None,
) -> list[GeneRecurrence]:
    """Fetch every gene with an SV record in ``study_id``, with its
    cohort-wide distinct-patient recurrence count, in one API call."""
    records = cbioportal_api.fetch_structural_variant_genes(
        [study_id], base_url=base_url, session=session
    )
    recurrence: list[GeneRecurrence] = []
    for record in records:
        symbol = record.get("hugoGeneSymbol")
        if not symbol:
            continue
        recurrence.append(
            GeneRecurrence(
                hugo_gene_symbol=symbol,
                entrez_gene_id=record.get("entrezGeneId"),
                distinct_patient_count=int(record.get("numberOfAlteredCases") or 0),
                total_sv_count=int(record.get("totalCount") or 0),
            )
        )
    return recurrence


@dataclass
class RecurrenceGateResult:
    """The result of applying a minimum-recurrence gate to a cohort's genes.

    ``total_genes`` and ``filtered_out_genes`` are always populated
    alongside ``passing_genes`` -- the pre-gate universe of genes is never
    dropped from this result, only from the downstream per-gene analysis.
    """

    study_id: str
    min_distinct_patients: int
    total_genes: int
    passing_genes: list[GeneRecurrence] = field(default_factory=list)
    filtered_out_genes: list[GeneRecurrence] = field(default_factory=list)

    @property
    def passing_count(self) -> int:
        return len(self.passing_genes)

    @property
    def filtered_out_count(self) -> int:
        return len(self.filtered_out_genes)


def gate_genes_by_recurrence(
    recurrence: list[GeneRecurrence],
    *,
    min_distinct_patients: int = DEFAULT_MIN_DISTINCT_PATIENTS,
    study_id: str = "",
) -> RecurrenceGateResult:
    """Split ``recurrence`` into genes meeting/not meeting the patient-count gate.

    Passing genes are sorted by descending recurrence (ties broken by gene
    symbol) so downstream consumers see the most-recurrent candidates
    first. Every gene in ``recurrence`` ends up in exactly one of
    ``passing_genes``/``filtered_out_genes`` -- ``total_genes`` is the
    count before this split was applied, so a caller can always report
    "N of M genes passed the >=k-patient gate" rather than only the
    post-filter count.
    """
    if min_distinct_patients < 0:
        raise ValueError(
            f"min_distinct_patients must be non-negative; got {min_distinct_patients!r}"
        )

    passing = [
        gene for gene in recurrence if gene.distinct_patient_count >= min_distinct_patients
    ]
    filtered_out = [
        gene for gene in recurrence if gene.distinct_patient_count < min_distinct_patients
    ]
    passing.sort(key=lambda gene: (-gene.distinct_patient_count, gene.hugo_gene_symbol))
    return RecurrenceGateResult(
        study_id=study_id,
        min_distinct_patients=min_distinct_patients,
        total_genes=len(recurrence),
        passing_genes=passing,
        filtered_out_genes=filtered_out,
    )

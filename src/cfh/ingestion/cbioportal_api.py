"""Client for the cBioPortal structural-variant REST API.

Gene selection (which Entrez ids to fetch) is always caller-supplied --
typically from a ``GeneConfig``'s ``entrez_gene_id`` -- so this module
stays generic across genes. The real-network call in
:func:`fetch_structural_variants` is only ever exercised by tests marked
``@pytest.mark.network`` (excluded from the default ``pytest`` run);
everything else in this module is plain, mockable request-building logic.
"""

from __future__ import annotations

import time
from typing import Iterable

import pandas as pd
import requests

from cfh.ingestion.sv_parser import OUTPUT_COLUMNS

DEFAULT_BASE_URL = "https://www.cbioportal.org/api"
DEFAULT_STUDY_ID = "msk_impact_50k_2026"
DEFAULT_SV_MOLECULAR_PROFILE_ID = "msk_impact_50k_2026_structural_variants"
_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}

_API_TO_NORMALIZED_COLUMNS = {
    "sampleId": "Sample_Id",
    "site1HugoSymbol": "Site1_Hugo_Symbol",
    "site1Chromosome": "Site1_Chromosome",
    "site1Position": "Site1_Position",
    "site2HugoSymbol": "Site2_Hugo_Symbol",
    "site2Chromosome": "Site2_Chromosome",
    "site2Position": "Site2_Position",
    "site2EffectOnFrame": "Site2_Effect_On_Frame",
    "tumorSplitReadCount": "Tumor_Split_Read_Count",
    "tumorPairedEndReadCount": "Tumor_Paired_End_Read_Count",
    "svStatus": "SV_Status",
    "ncbiBuild": "NCBI_Build",
    "connectionType": "Connection_Type",
    "breakpointType": "Breakpoint_Type",
    "annotation": "Annotation",
    "eventInfo": "Event_Info",
}


def fetch_structural_variants(
    entrez_gene_ids: Iterable[int],
    molecular_profile_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 3,
    backoff_seconds: float = 0.5,
) -> list[dict]:
    """POST to ``/structural-variant/fetch`` and return the parsed JSON body.

    ``molecular_profile_ids`` is required and has no default: which cohort's
    SV profile to query is always caller-supplied (e.g. from ingestion
    config), never silently defaulted to a specific study like MSK-IMPACT.
    ``DEFAULT_SV_MOLECULAR_PROFILE_ID`` remains available for callers that
    do want the MSK-IMPACT 50k profile, but it's opt-in, not automatic.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/structural-variant/fetch"
    body = {
        "entrezGeneIds": list(entrez_gene_ids),
        "molecularProfileIds": list(molecular_profile_ids),
    }
    attempt = 0
    while True:
        response = session.post(url, json=body, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(backoff_seconds * (2**attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def fetch_structural_variant_genes(
    study_ids: Iterable[str],
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
    max_retries: int = 3,
    backoff_seconds: float = 0.5,
) -> list[dict]:
    """POST to ``/structuralvariant-genes/fetch`` for cohort-wide SV gene recurrence.

    Returns one ``AlterationCountByGene``-shaped dict per gene that has at
    least one structural-variant record anywhere in the requested
    studies -- every such gene in one call, each carrying its own
    ``numberOfAlteredCases`` (distinct-patient count) and ``totalCount``
    (raw SV record count). This is the whole-cohort recurrence signal used
    to gate a genome-wide scan down to recurrently-altered genes, as
    opposed to :func:`fetch_structural_variants`, which fetches per-event
    SV records for an already-selected set of Entrez gene ids.
    """
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/structuralvariant-genes/fetch"
    body = {"studyIds": list(study_ids)}
    attempt = 0
    while True:
        response = session.post(url, json=body, timeout=timeout)
        if response.status_code not in _RETRYABLE_STATUS_CODES or attempt >= max_retries:
            break
        time.sleep(backoff_seconds * (2**attempt))
        attempt += 1
    response.raise_for_status()
    return response.json()


def structural_variants_to_dataframe(calls: Iterable[dict]) -> pd.DataFrame:
    """Adapt cBioPortal camelCase API objects to the production SV schema."""
    records = []
    for row_number, call in enumerate(calls, start=1):
        record = {
            destination: call.get(source)
            for source, destination in _API_TO_NORMALIZED_COLUMNS.items()
        }
        record["Extra_fields"] = {
            key: value for key, value in call.items() if key not in _API_TO_NORMALIZED_COLUMNS
        }
        record["Source_row_number"] = row_number
        record["Parse_warnings"] = None
        records.append(record)
    return pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)

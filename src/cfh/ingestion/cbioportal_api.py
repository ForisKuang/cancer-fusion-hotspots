"""Client for the cBioPortal structural-variant REST API.

Gene selection (which Entrez ids to fetch) is always caller-supplied --
typically from a ``GeneConfig``'s ``entrez_gene_id`` -- so this module
stays generic across genes. The real-network call in
:func:`fetch_structural_variants` is only ever exercised by tests marked
``@pytest.mark.network`` (excluded from the default ``pytest`` run);
everything else in this module is plain, mockable request-building logic.
"""

from __future__ import annotations

from typing import Iterable

import requests

DEFAULT_BASE_URL = "https://www.cbioportal.org/api"
DEFAULT_STUDY_ID = "msk_impact_50k_2026"
DEFAULT_SV_MOLECULAR_PROFILE_ID = "msk_impact_50k_2026_structural_variants"


def fetch_structural_variants(
    entrez_gene_ids: Iterable[int],
    molecular_profile_ids: Iterable[str] = (DEFAULT_SV_MOLECULAR_PROFILE_ID,),
    *,
    base_url: str = DEFAULT_BASE_URL,
    session: "requests.Session | None" = None,
    timeout: float = 30,
) -> list[dict]:
    """POST to ``/structural-variant/fetch`` and return the parsed JSON body."""
    session = session or requests.Session()
    url = f"{base_url.rstrip('/')}/structural-variant/fetch"
    body = {
        "entrezGeneIds": list(entrez_gene_ids),
        "molecularProfileIds": list(molecular_profile_ids),
    }
    response = session.post(url, json=body, timeout=timeout)
    response.raise_for_status()
    return response.json()

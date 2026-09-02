from unittest.mock import MagicMock

import pytest

from cfh.ingestion import cbioportal_api


def test_fetch_structural_variants_posts_expected_body_without_real_network():
    mock_session = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = [{"sampleId": "SAMPLE-001"}]
    mock_response.raise_for_status.return_value = None
    mock_session.post.return_value = mock_response

    result = cbioportal_api.fetch_structural_variants(
        [cbioportal_api.BRAF_ENTREZ_GENE_ID],
        [cbioportal_api.DEFAULT_SV_MOLECULAR_PROFILE_ID],
        session=mock_session,
    )

    assert result == [{"sampleId": "SAMPLE-001"}]
    mock_session.post.assert_called_once()
    _, kwargs = mock_session.post.call_args
    called_url = mock_session.post.call_args.args[0]
    assert called_url == f"{cbioportal_api.DEFAULT_BASE_URL}/structural-variant/fetch"
    assert kwargs["json"] == {
        "entrezGeneIds": [673],
        "molecularProfileIds": [cbioportal_api.DEFAULT_SV_MOLECULAR_PROFILE_ID],
    }


def test_fetch_braf_structural_variants_uses_defaults_without_real_network():
    mock_session = MagicMock()
    mock_response = MagicMock()
    mock_response.json.return_value = []
    mock_session.post.return_value = mock_response

    cbioportal_api.fetch_braf_structural_variants(session=mock_session)

    _, kwargs = mock_session.post.call_args
    assert kwargs["json"]["entrezGeneIds"] == [673]
    assert kwargs["json"]["molecularProfileIds"] == [
        "msk_impact_50k_2026_structural_variants"
    ]


@pytest.mark.network
def test_fetch_structural_variants_real_network_call():
    """Excluded from default `pytest -m "not network"` runs."""
    result = cbioportal_api.fetch_braf_structural_variants()
    assert isinstance(result, list)

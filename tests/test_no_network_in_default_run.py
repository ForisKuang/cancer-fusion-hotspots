"""Positive-proof test: the repo-wide network guard (see
``conftest.py::_block_real_network_calls``, an autouse fixture applied to
every test not marked ``@pytest.mark.network``) actually blocks a real
outbound call, rather than the mocked call underneath it just happening to
never touch the network.
"""

from unittest.mock import MagicMock

from cfh.ingestion import cbioportal_api


def test_cbioportal_client_mocked_call_makes_no_real_socket():
    mock_session = MagicMock()
    mock_session.post.return_value.json.return_value = []
    cbioportal_api.fetch_structural_variants([673], session=mock_session)
    # No AssertionError raised above means no real socket was opened
    # (the autouse conftest fixture would have raised one otherwise).

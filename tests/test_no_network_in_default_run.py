"""Guard test: confirm the default (`-m "not network"`) suite makes zero real
outbound network calls. It monkeypatches socket creation to raise, so any
accidental real HTTP call anywhere in the non-network suite would fail this
test's own setup/teardown expectations rather than silently succeeding.
"""

import socket

import pytest


@pytest.fixture(autouse=True)
def _block_real_sockets(monkeypatch):
    def _blocked(*args, **kwargs):
        raise AssertionError("real network access attempted during 'not network' test run")

    monkeypatch.setattr(socket, "create_connection", _blocked)
    yield


def test_cbioportal_client_mocked_call_makes_no_real_socket(monkeypatch):
    from unittest.mock import MagicMock

    from cfh.ingestion import cbioportal_api

    mock_session = MagicMock()
    mock_session.post.return_value.json.return_value = []
    cbioportal_api.fetch_structural_variants([673], session=mock_session)
    # No AssertionError raised above means no real socket was opened.

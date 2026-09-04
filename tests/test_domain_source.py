import json
from unittest.mock import MagicMock

from cfh.mapping.domain_source import UniProtDomainSource

# Fixture frozen from a real response:
#   GET https://rest.uniprot.org/uniprotkb/P15056.json
# fetched 2026-09-02. See tests/fixtures/uniprot/P15056.json.


def _load_fixture(uniprot_fixture_path):
    return json.loads(uniprot_fixture_path.read_text())


def test_parse_returns_kinase_domain_with_fixtures_own_coordinates(uniprot_fixture_path):
    payload = _load_fixture(uniprot_fixture_path)
    domains = UniProtDomainSource().parse(payload)

    kinase_domains = [d for d in domains if "protein kinase" in d.name.lower()]
    assert len(kinase_domains) == 1
    kinase = kinase_domains[0]

    # Assert against the fixture's own parsed values, not hardcoded literature numbers.
    expected_feature = next(
        f
        for f in payload["features"]
        if f.get("type") == "Domain" and "kinase" in (f.get("description") or "").lower()
    )
    assert kinase.start_aa == expected_feature["location"]["start"]["value"]
    assert kinase.end_aa == expected_feature["location"]["end"]["value"]


def test_fetch_caches_after_first_call_one_http_call_total(uniprot_fixture_path):
    mock_session = MagicMock()
    mock_session.get.return_value.json.return_value = _load_fixture(uniprot_fixture_path)
    mock_session.get.return_value.raise_for_status.return_value = None

    source = UniProtDomainSource(session=mock_session)
    first = source.fetch("P15056")
    second = source.fetch("P15056")

    assert first == second
    mock_session.get.assert_called_once()


def test_fetch_uses_on_disk_cache_without_repopulating_http(tmp_path, uniprot_fixture_path):
    cache_dir = tmp_path / "domain_cache"
    cache_dir.mkdir()
    (cache_dir / "P15056.json").write_text(uniprot_fixture_path.read_text())

    mock_session = MagicMock()
    source = UniProtDomainSource(session=mock_session, cache_dir=cache_dir)

    domains = source.fetch("P15056")

    assert any("kinase" in d.name.lower() for d in domains)
    mock_session.get.assert_not_called()

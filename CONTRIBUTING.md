# Contributing to Cancer Fusion Hotspots

Thanks for your interest in contributing.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Running checks locally

```bash
ruff check .
pytest -m "not network"
```

Tests marked `@pytest.mark.network` make real calls to external services
(cBioPortal, Ensembl, UniProt, InterPro) and are excluded from the default
test run. CI never runs them; run them locally only when you need to
validate against a live upstream API.

## Adding support for a new gene

Gene biology (canonical transcript, protein accession, key domains,
autoinhibitory domains) belongs in a new YAML file under
`src/cfh/genes/configs/`, loaded via `GeneConfig`. Generic
ingestion/normalization/mapping code must never hardcode a specific gene
symbol, transcript ID, or protein accession — those always come from the
`GeneConfig` passed in by the caller.

## Adding a hotspot-detection algorithm

Implement the `Algorithm` interface in `src/cfh/algorithms/base.py` and
register it via `src/cfh/algorithms/registry.py` so it can be discovered by
name.

## Pull requests

- Keep changes focused and covered by tests.
- Run `ruff check .` and `pytest -m "not network"` before opening a PR.
- Describe what changed and how it was verified.

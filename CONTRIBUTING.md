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
ruff format --check .
mypy src
pytest -m "not network"
```

The same lint and formatting checks can run automatically before each commit:

```bash
pre-commit install
pre-commit run --all-files
```

Tests marked `@pytest.mark.network` make real calls to external services
(cBioPortal, Ensembl, UniProt, InterPro) and are excluded from the default
test run. A weekly and manually dispatchable workflow runs them against the
live services. To run them locally, opt in explicitly:

```bash
CFH_RUN_NETWORK_TESTS=1 pytest -m network
```

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
- Run `pre-commit run --all-files`, `mypy src`, and
  `pytest -m "not network"` before opening a PR.
- Describe what changed and how it was verified.

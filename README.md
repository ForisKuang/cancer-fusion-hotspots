# Cancer Fusion Hotspots (cfh)

`cfh` detects gene-fusion breakpoint and domain-retention hotspots from
structural-variant calls in cancer sequencing cohorts. Milestone 1 ("Stage A:
core pipeline") targets **BRAF** fusions in the public MSK-IMPACT 50k cohort
on cBioPortal, using a design that generalizes to other genes without
code changes.

> This repository is currently a placeholder staged for eventual transfer to
> `cBioPortal/cancer-fusion-hotspots`.

## What this is

Gene fusions can retain, disrupt, or lose functional protein domains
depending on exactly where the breakpoint falls. `cfh` builds a reproducible
pipeline that:

1. **Ingests** structural-variant (SV) and clinical data from cBioPortal,
   either from a downloaded study archive or the cBioPortal REST API.
2. **Normalizes** raw SV rows into a typed `FusionEvent` model, classifying
   event type (fusion / deletion / inversion / translocation), reading
   frame, 5'/3' orientation, and tumor-type/OncoTree provenance (joined
   from clinical data) without ever guessing when the source data is
   ambiguous.
3. **Maps** each event's breakpoint onto transcript exons and protein
   domains (via RefSeq annotations, Genome Nexus's canonical-transcript
   and exon-coordinate data as the primary generic fallback, and an
   Ensembl REST / UniProt+InterPro cross-check), producing a
   `FusionFeature` describing which domains are retained, disrupted, or
   lost.
4. **Runs pluggable hotspot-detection algorithms** against the normalized
   events and features, each producing a structured `AlgorithmResult`.

Gene-specific biology (canonical transcript, protein accession, domain
boundaries) lives in a per-gene YAML config (see
`src/cfh/genes/configs/`), never hardcoded into the generic
ingestion/normalization/mapping code, so adding a new gene is a config
change, not a code change.

## Install

Requires Python 3.10+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Quickstart

```python
from cfh.genes.registry import load_gene_config
from cfh.ingestion.archive_reader import load_sv_dataframe
from cfh.ingestion.clinical_parser import parse_clinical_sample
from cfh.normalization.event_normalizer import normalize

gene_config = load_gene_config("braf")

raw_sv = load_sv_dataframe("path/to/msk_impact_50k_study")
clinical = parse_clinical_sample("path/to/msk_impact_50k_study/data_clinical_sample.txt")

events = normalize(raw_sv, clinical, cohort="msk_impact_50k_2026")
print(f"Normalized {len(events)} fusion events for {gene_config.gene_symbol}")
```

Run the test suite (excluding tests that hit real external services):

```bash
pytest -m "not network"
```

Run the live cBioPortal/Genome Nexus benchmark. Each invocation creates a
timestamped directory under `runs/` with a manifest, event-level TSV, structured
JSON, Markdown report, discrepancy table, and SVG visualizations:

```bash
cfh real-benchmark BRAF msk_impact_50k_2026
```

The study ID is an argument to the shared pipeline. For example, the original
MSK-IMPACT publication cohort can be run with:

```bash
cfh real-benchmark BRAF msk_impact_2017
```

To run the complete registered-algorithm orchestrator for a configured gene and
study, use:

```bash
cfh analyze BRAF msk_impact_50k_2026
```

RET uses the same command and live ingestion/mapping path:

```bash
cfh analyze RET msk_impact_50k_2026
```

The TCGA PanCancer Atlas studies expose structural variants under the same
`<study_id>_structural_variants` profile convention. The BRAF holdout uses the
thyroid carcinoma study, which contains BRAF fusion calls:

```bash
cfh real-benchmark BRAF thca_tcga_pan_can_atlas_2018
```

All 32 PanCancer Atlas study IDs with `data_sv.txt`/`meta_sv.txt` in the public
cBioPortal Datahub are: `acc`, `blca`, `brca`, `cesc`, `chol`, `coadread`,
`dlbc`, `esca`, `gbm`, `hnsc`, `kich`, `kirc`, `kirp`, `laml`, `lgg`, `lihc`,
`luad`, `lusc`, `meso`, `ov`, `paad`, `pcpg`, `prad`, `sarc`, `skcm`, `stad`,
`tgct`, `thca`, `thym`, `ucec`, `ucs`, and `uvm`, each followed by
`_tcga_pan_can_atlas_2018`.

These commands make unauthenticated requests to both public services.

## Repository layout

```
src/cfh/model/         FusionEvent, FusionFeature, AlgorithmResult schemas
src/cfh/genes/         Per-gene configuration registry (e.g. braf.yaml)
src/cfh/algorithms/    Hotspot-detection algorithm plugin interface + registry
src/cfh/ingestion/     cBioPortal archive/API ingestion
src/cfh/normalization/ Raw SV rows -> FusionEvent
src/cfh/mapping/       Transcript/exon/domain mapping -> FusionFeature
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT — see [LICENSE](LICENSE).

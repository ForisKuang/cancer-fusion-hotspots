# BRAF real-data fusion benchmark: msk_impact_2017

Retrieved from public cBioPortal and Genome Nexus on 2026-09-03.

## Results

- Structural variants returned for BRAF: 48
- Protein-fusion records found: 35
- Protein-fusion records mapped: 35
- Malformed/unmappable fusion records skipped: 0
- In-frame: 33/35 (94.3%)
- PF07714 (458-712 aa) retained: 31/35 (88.6%)
- In-frame and PF07714-retained: 31/33
- Fisher exact test (one-sided): odds ratio inf, p=0.010084
- Breakpoint-permutation empirical p-value: 0.342657
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[31, 0], [2, 2]]`

## Method

The cBioPortal `msk_impact_2017_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Records explicitly annotated as `Protein Fusion` were adapted to the production SV schema and normalized; the raw `site2EffectOnFrame=NA` values were therefore resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

BRAF genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

## Partners

AGAP3 (2), AGK (3), CCDC6 (1), CDK5RAP2 (2), CUL1 (1), FAM131B (1), GIPC2 (1), KIAA1549 (4), MKRN1 (3), OSBPL9 (1), PARP12 (1), PHTF2 (1), PJA2 (1), PRKAR1B (1), PRKAR2B (1), RBM33 (1), SCRIB (1), SND1 (8), ZNF207 (1)

## Reference comparison

| Metric | PMC5461196 | This run |
|---|---:|---:|
| In-frame | 100.0% | 94.3% |
| Domain retained | 100.0% | 88.6% |

## Interpretation

These values describe the live study named above.

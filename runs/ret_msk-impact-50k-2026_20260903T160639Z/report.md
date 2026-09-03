# RET real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-03.

## Results

- Structural variants returned for RET: 230
- Protein-fusion records found: 194
- Protein-fusion records mapped: 194
- Malformed/unmappable fusion records skipped: 0
- In-frame: 146/194 (75.3%)
- PF07714 (724-1005 aa) retained: 144/194 (74.2%)
- In-frame and PF07714-retained: 119/146
- Fisher exact test (one-sided): odds ratio 4.05481, p=9.84123e-05
- Breakpoint-permutation empirical p-value: 0.00699301
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[119, 25], [27, 23]]`

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Fusion-annotated records were adapted to the production SV schema and normalized; when `site2EffectOnFrame=NA`, frame status was resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

RET genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

## Partners

CCDC186 (1), CCDC6 (51), CLIP1 (1), CTNNA1 (1), CTNNA3 (1), CUBN (1), DOCK1 (1), EML4 (1), ERC1 (1), ERCC6 (3), FBXL7 (2), GEMIN5 (1), GRIPAP1 (1), KCNMA1 (1), KIAA1217 (1), KIF13A (1), KIF5B (87), LIN52 (1), NCOA4 (16), PDCD10 (1), RASSF4 (1), RBPMS (1), RELCH (1), RUFY1 (2), RUFY2 (1), RUFY3 (2), SHROOM3 (1), SLC4A4 (1), SPECC1L (4), TFG (1), TIMM23B (3), TRIM24 (1), UXS1 (1)

## Reference comparison

| Metric | PMCID: PMC6430196 (KIF5B-RET) | This run |
|---|---:|---:|
| In-frame | 100.0% | 75.3% |
| Domain retained | 100.0% | 74.2% |

## Interpretation

These values describe the live study named above.

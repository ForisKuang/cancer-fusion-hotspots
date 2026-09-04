# BRAF real-data fusion benchmark: msk_impact_50k_2026

Retrieved from public cBioPortal and Genome Nexus on 2026-09-03.

## Results

- Structural variants returned for BRAF: 251
- Protein-fusion records found: 174
- Protein-fusion records mapped: 174
- Malformed/unmappable fusion records skipped: 0
- In-frame: 151/174 (86.8%)
- PF07714 (458-712 aa) retained: 161/174 (92.5%)
- In-frame and PF07714-retained: 142/151
- Fisher exact test (one-sided): odds ratio 3.32164, p=0.0737799
- Breakpoint-permutation empirical p-value: 0.0719281
- Contingency table `[[retained/in-frame, retained/other], [not-retained/in-frame, not-retained/other]]`: `[[142, 19], [9, 4]]`

## Method

The cBioPortal `msk_impact_50k_2026_structural_variants` structural-variant profile was queried by the configured Entrez gene ID. Records explicitly annotated as `Protein Fusion` were adapted to the production SV schema and normalized; the raw `site2EffectOnFrame=NA` values were therefore resolved from `Event_Info`, not copied into `FusionEvent.Frame_status`.

BRAF genomic breakpoints were mapped against the Genome Nexus canonical transcript, and retention was classified against its returned PF07714 coordinates. Counts are event-level with no patient deduplication. The Fisher comparison's `other` column combines out-of-frame and unknown-frame events, as pre-specified by the domain-retention algorithm.

## Partners

ABCC1 (1), ABCC2 (1), AGAP3 (6), AGK (14), AKAP9 (2), ATF7 (1), CAPZA2 (1), CARMIL1 (2), CCDC6 (1), CDK5RAP2 (4), CLEC2L (1), CLIP1 (2), CUL1 (5), DBI (1), DOCK4 (1), EPB41 (1), ERC1 (1), EXOC4 (2), FAM131B (3), GIPC2 (1), GTF2I (1), GTF2IRD1 (3), ITGB5 (1), KCND2 (1), KIAA1549 (43), KLHL12 (1), LRBA (1), MAD1L1 (2), MINDY4 (1), MKLN1 (1), MKRN1 (10), NEO1 (1), NFIA (1), NRF1 (2), OSBPL7 (1), OSBPL9 (1), PAK2 (1), PHTF2 (1), PJA2 (1), PLPP1 (1), PRIM2 (1), PRKAR1B (1), PRKAR2B (1), PTPRZ1 (1), PWWP2A (1), RBM33 (1), RBMS3 (1), RGS3 (1), RRBP1 (1), SCRIB (1), SETBP1 (1), SND1 (17), SNX8 (1), SPOCK1 (1), TBXAS1 (1), TMPRSS2 (3), TRA2B (1), TRIM24 (8), TRIM33 (1), UBE2H (1), WEE2 (1), ZC3HAV1 (1), ZNF207 (1)

## Reference comparison

| Metric | PMC5461196 | This run |
|---|---:|---:|
| In-frame | 100.0% | 86.8% |
| Domain retained | 100.0% | 92.5% |

## Interpretation

This does **not** reproduce the Zehir et al. (PMC5461196) report of 33/33 BRAF fusions being in-frame with the kinase domain retained: this live successor cohort has 151/174 in-frame and 142/151 in-frame fusions retaining PF07714.

`msk_impact_50k_2026` is a newer successor cohort, not the paper's original `msk_impact_2017` cohort. This is therefore replication in a related cohort, not a reanalysis of the paper's original 33 cases.

# Genome-wide fusion-hotspot cohort scan: msk_impact_50k_2026

- Total genes with any structural-variant record in the cohort: 3919
- Genes passing the >= 5-distinct-patient recurrence gate: 544
- Curated gene configs used: 2
- Auto-generated gene configs used: 522
- Genes gated in but unresolvable (no Genome Nexus canonical transcript): 20
- FDR-significant genes (q < 0.05) after Benjamini-Hochberg correction across all 544 scanned genes: 2

## Warnings

- 20 gated gene(s) had no resolvable canonical transcript in Genome Nexus and were skipped: AGAP3, CDH3, CDKN2B-AS1, CREB3L1, ERF, H3C2, H3C6, KMT2B, LINC00114, MAP3K14, MEF2B, NBPF20, PLCD3, RAB11FIP4, RECQL4, SEPTIN14, TCF3, TRAF2, TRAP1, ZFTA

## Scanned genes (sorted by significance)

| gene_symbol | config_source | status | distinct_patient_count | total_sv_count | n_events_analyzed | in_frame_percent | domain_retention_percent | fisher_p_value | permutation_p_value | min_fdr_adjusted_q_value | fdr_significant | top_composite_score | top_composite_partner_gene | error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ETV6 | auto | ok | 213 | 214 | 90 | 71.11 | 75.56 | 5.123e-06 | 0.009901 | 0.004437 | yes | 0.4164 | NTRK3 |  |
| RET | curated | ok | 219 | 230 | 194 | 75.26 | 74.23 | 9.841e-05 | 0.009901 | 0.04261 | yes | 0.4045 | KIF5B |  |
| TMPRSS2 | auto | ok | 1088 | 1132 | 867 | 29.99 | 43.37 | 0.4193 | 0.08192 | 0.1453 | no | 0.5149 | ERG |  |
| ERG | auto | ok | 863 | 863 | 788 | 30.2 | 96.07 | 0.8462 | 0.4653 | 0.1453 | no | 0.4294 | TMPRSS2 |  |
| EGFR | auto | ok | 466 | 559 | 55 | 50.91 | 63.64 | 0.01893 | 0.6337 | 0.1453 | no | 0.2369 | RAD51 |  |
| ALK | auto | ok | 321 | 322 | 272 | 83.09 | 95.96 | 0.001917 | 0.009901 | 0.1453 | no | 0.5133 | EML4 |  |
| TP53 | auto | ok | 258 | 260 | 47 | 8.511 | 31.91 | 0.8341 | 1 | 0.1453 | no | 0.2508 | EIF5 |  |
| BRAF | curated | ok | 247 | 251 | 179 | 84.36 | 91.06 | 0.01337 | 0.04995 | 0.1453 | no | 0.2681 | KIAA1549 |  |
| EML4 | auto | ok | 224 | 225 | 223 | 87.89 | 29.15 | 0.993 | 0.6931 | 0.1453 | no | 0.558 | ALK |  |
| ROS1 | auto | ok | 204 | 207 | 122 | 59.02 | 86.07 | 0.0728 | 0.009901 | 0.1453 | no | 0.3779 | CD74 |  |
| FGFR3 | auto | ok | 194 | 195 | 152 | 48.68 | 39.47 | 0.4872 | 0.9901 | 0.1453 | no | 0.4396 | TACC3 |  |
| FGFR2 | auto | ok | 188 | 195 | 136 | 79.41 | 86.03 | 0.005225 | 0.009901 | 0.1453 | no | 0.3665 | BICC1 |  |
| RB1 | auto | ok | 147 | 147 | 19 | 0 | 47.37 |  |  | 0.1453 | no | 0.3684 | PHF11 |  |
| TACC3 | auto | ok | 138 | 139 | 132 | 50 | 26.52 | 0.8817 | 1 | 0.1453 | no | 0.4766 | FGFR3 |  |
| CREBBP | auto | ok | 121 | 122 | 46 | 13.04 | 60.87 | 0.04022 | 0.03596 | 0.1453 | no | 0.3224 | TRAP1 |  |
| WT1 | auto | ok | 119 | 119 | 104 | 89.42 | 6.731 | 0.9999 | 0.9901 | 0.1453 | no | 0.6177 | EWSR1 |  |
| NTRK1 | auto | ok | 113 | 118 | 78 | 41.03 | 71.79 | 0.003635 | 0.009901 | 0.1453 | no | 0.3245 | LMNA |  |
| NAB2 | auto | ok | 111 | 111 | 72 | 34.72 | 68.06 | 0.1126 | 0.009901 | 0.1453 | no | 0.4883 | STAT6 |  |
| CDK12 | auto | ok | 99 | 101 | 43 | 16.28 | 18.6 | 0.81 | 0.9802 | 0.1453 | no | 0.2345 | ADCY9 |  |
| TERT | auto | ok | 92 | 92 | 19 | 0 | 89.47 |  |  | 0.1453 | no | 0.3988 | SLC12A7 |  |
| CD74 | auto | ok | 74 | 74 | 67 | 58.21 | 88.06 | 0.2987 | 0.009901 | 0.1453 | no | 0.549 | ROS1 |  |
| MET | auto | ok | 73 | 75 | 33 | 51.52 | 60.61 | 0.08482 | 0.009901 | 0.1453 | no | 0.3118 | CNTNAP2 |  |
| NTRK3 | auto | ok | 71 | 71 | 58 | 94.83 | 68.97 | 1 | 1 | 0.1453 | no | 0.5481 | ETV6 |  |
| STAT6 | auto | ok | 60 | 60 | 53 | 43.4 | 11.32 | 0.8314 | 1 | 0.1453 | no | 0.4286 | NAB2 |  |
| DNAJB1 | auto | ok | 54 | 54 | 39 | 64.1 | 30.77 | 0.01729 | 1 | 0.1453 | no | 0.4259 | PRKACA |  |
| CCDC6 | auto | ok | 53 | 53 | 53 | 92.45 | 43.4 | 0.7883 | 0.009901 | 0.1453 | no | 0.484 | RET |  |
| TFE3 | auto | ok | 45 | 46 | 37 | 56.76 | 75.68 | 0.1714 | 0.6139 | 0.1453 | no | 0.2946 | ASPSCR1 |  |
| KIAA1549 | auto | ok | 44 | 44 | 44 | 90.91 | 50 | 0.9461 | 0.009901 | 0.1453 | no | 0.6367 | BRAF |  |
| POLD1 | auto | ok | 36 | 37 | 16 | 12.5 | 37.5 | 0.1923 | 0.4356 | 0.1453 | no | 0.2807 | MYH14 |  |
| RBM10 | auto | ok | 36 | 36 | 7 | 42.86 | 71.43 | 0.2857 | 0.009901 | 0.1453 | no | 0.3471 | PHIP |  |
| PRKACA | auto | ok | 35 | 35 | 35 | 68.57 | 97.14 | 0.3143 | 0.009901 | 0.1453 | no | 0.6959 | DNAJB1 |  |
| ATF1 | auto | ok | 30 | 30 | 30 | 70 | 90 | 0.6724 | 0.009901 | 0.1453 | no | 0.4859 | EWSR1 |  |
| ETV1 | auto | ok | 30 | 30 | 21 | 52.38 | 4.762 | 0.5238 | 0.009901 | 0.1453 | no | 0.516 | TMPRSS2 |  |
| BICC1 | auto | ok | 29 | 29 | 28 | 96.43 | 35.71 | 0.6429 | 1 | 0.1453 | no | 0.4286 | FGFR2 |  |
| SLX4 | auto | ok | 28 | 29 | 19 | 26.32 | 68.42 | 0.5675 | 0.1287 | 0.1453 | no | 0.2896 | CREBBP |  |
| LMNA | auto | ok | 17 | 17 | 13 | 46.15 | 38.46 | 0.4126 | 0.06194 | 0.1453 | no | 0.4387 | NTRK1 |  |
| NRG1 | auto | ok | 16 | 16 | 15 | 53.33 | 60 | 0.9161 | 0.9802 | 0.1453 | no | 0.4003 | CD74 |  |
| TTC28 | auto | ok | 15 | 15 | 12 | 16.67 | 58.33 | 0.3818 | 0.009901 | 0.1453 | no | 0.4144 | CHEK2 |  |
| AGK | auto | ok | 14 | 14 | 14 | 92.86 | 14.29 | 0.8571 | 0.009901 | 0.1453 | no | 0.5002 | BRAF |  |
| TPM3 | auto | ok | 14 | 14 | 14 | 78.57 | 64.29 | 1 | 0.009901 | 0.1453 | no | 0.5289 | NTRK1 |  |
| EZR | auto | ok | 13 | 13 | 13 | 61.54 | 30.77 | 0.4895 | 0.009901 | 0.1453 | no | 0.5002 | ROS1 |  |
| TRAPPC9 | auto | ok | 13 | 13 | 13 | 23.08 | 30.77 | 0.7063 | 0.03796 | 0.1453 | no | 0.5743 | AGO2 |  |
| GOPC | auto | ok | 6 | 6 | 4 | 75 | 75 | 0.25 | 0.009901 | 0.1453 | no | 0.5002 | ROS1 |  |
| EMID1 | auto | ok | 5 | 5 | 4 | 75 | 25 | 0.75 | 0.009901 | 0.1453 | no | 0.5002 | EWSR1 |  |
| AKT2 | auto | ok | 22 | 22 | 10 | 10 | 40 | 1 | 1 | 0.2372 | no | 0.2571 | PRR5L |  |
| NCOR1 | auto | ok | 80 | 81 | 19 | 0 | 36.84 |  |  | 0.2437 | no | 0.3684 | CUX1 |  |
| NCOA3 | auto | ok | 33 | 34 | 13 | 23.08 | 46.15 | 0.5 | 0.0198 | 0.2437 | no | 0.333 | EYA2 |  |
| SRP19 | auto | ok | 15 | 15 | 8 | 0 | 62.5 |  |  | 0.2437 | no | 0.6667 | APC |  |
| FOXO1 | auto | ok | 14 | 14 | 4 | 50 | 25 | 0.5 | 0.0198 | 0.2437 | no | 0.4507 | PAX3 |  |
| STRN | auto | ok | 8 | 8 | 8 | 25 | 50 | 0.2143 | 0.02498 | 0.2437 | no | 0.4322 | ALK |  |
| CREB1 | auto | ok | 7 | 7 | 7 | 57.14 | 71.43 | 0.7143 | 0.0198 | 0.2437 | no | 0.4894 | EWSR1 |  |
| INPPL1 | auto | ok | 39 | 39 | 14 | 28.57 | 42.86 | 0.02098 | 0.1089 | 0.2523 | no | 0.2793 | ANK2 |  |
| RNF43 | auto | ok | 48 | 48 | 21 | 28.57 | 42.86 | 0.9344 | 1 | 0.2607 | no | 0.2347 | EFCAB5 |  |
| DNMT1 | auto | ok | 80 | 82 | 28 | 17.86 | 28.57 | 0.8835 | 0.4653 | 0.323 | no | 0.275 | RAVER1 |  |
| FLI1 | auto | ok | 121 | 121 | 118 | 61.86 | 4.237 | 0.9935 | 0.6634 | 0.3882 | no | 0.6142 | EWSR1 |  |
| NOTCH1 | auto | ok | 115 | 117 | 47 | 23.4 | 59.57 | 0.5602 | 0.6436 | 0.3893 | no | 0.3121 | SEC16A |  |
| AGO2 | auto | ok | 49 | 49 | 23 | 26.09 | 52.17 | 0.4171 | 0.03796 | 0.4009 | no | 0.4292 | TRAPPC9 |  |
| TOP1 | auto | ok | 31 | 31 | 10 | 0 | 50 |  |  | 0.4223 | no | 0.4267 | DHX35 |  |
| MAP3K1 | auto | ok | 42 | 42 | 10 | 10 | 30 | 1 | 1 | 0.458 | no | 0.2639 | PDE4D |  |
| EWSR1 | auto | ok | 415 | 418 | 349 | 65.33 | 7.45 | 0.9999 | 1 | 0.4916 | no | 0.2298 | FLI1 |  |
| TSC1 | auto | ok | 28 | 28 | 7 | 0 | 57.14 |  |  | 0.5249 | no | 0.4286 | TMC1 |  |
| FLT4 | auto | ok | 35 | 35 | 7 | 14.29 | 42.86 | 0.4286 | 0.505 | 0.5671 | no | 0.2861 | COL23A1 |  |
| NOTCH2 | auto | ok | 92 | 92 | 28 | 3.571 | 53.57 | 0.5769 | 0.1485 | 0.5799 | no | 0.2745 | CEP85 |  |
| SUFU | auto | ok | 17 | 17 | 8 | 25 | 25 | 0.4643 | 0.9208 | 0.5861 | no | 0.2691 | ETV6 |  |
| ATP1B2 | auto | ok | 7 | 7 | 6 | 0 | 66.67 |  |  | 0.5982 | no | 0.6667 | TP53 |  |
| FH | auto | ok | 29 | 30 | 16 | 6.25 | 12.5 | 0.1333 | 0.2871 | 0.6284 | no | 0.2604 | RYR2 |  |
| FBXL20 | auto | ok | 8 | 8 | 7 | 28.57 | 42.86 | 0.1429 | 0.07193 | 0.6489 | no | 0.4082 | CDK12 |  |
| ARID1B | auto | ok | 72 | 72 | 13 | 0 | 53.85 |  |  | 0.6554 | no | 0.5385 | GRIN2A |  |
| BAP1 | auto | ok | 67 | 70 | 18 | 16.67 | 61.11 | 0.674 | 0.1089 | 0.6554 | no | 0.291 | PBRM1 |  |
| EIF4E | auto | ok | 5 | 5 | 2 | 50 | 50 | 0.5 | 0.07692 | 0.6596 | no | 0.3234 | SHROOM3 |  |
| STAT5B | auto | ok | 27 | 27 | 13 | 15.38 | 69.23 | 0.9545 | 0.1188 | 0.6821 | no | 0.2803 | USP32 |  |
| STAT5A | auto | ok | 23 | 23 | 13 | 7.692 | 15.38 | 1 | 1 | 0.6821 | no | 0.2473 | EZH1 |  |
| BCAS3 | auto | ok | 18 | 19 | 16 | 6.25 | 43.75 | 1 | 1 | 0.6839 | no | 0.2411 | EZH1 |  |
| MTOR | auto | ok | 49 | 51 | 20 | 10 | 45 | 0.1895 | 0.08791 | 0.7115 | no | 0.2734 | CAPZB |  |
| NELL1 | auto | ok | 5 | 5 | 5 | 20 | 60 | 1 | 1 | 0.737 | no | 0.3857 | NOTCH1 |  |
| CDH1 | auto | ok | 58 | 58 | 9 | 11.11 | 55.56 | 0.7143 | 0.901 | 0.754 | no | 0.2635 | SMPD3 |  |
| ZNRF3 | auto | ok | 6 | 6 | 6 | 33.33 | 50 | 0.8 | 1 | 0.7708 | no | 0.4235 | EWSR1 |  |
| ATRX | auto | ok | 64 | 65 | 7 | 28.57 | 28.57 | 0.5238 | 0.1089 | 0.7731 | no | 0.3268 | MAGT1 |  |
| BRIP1 | auto | ok | 46 | 46 | 19 | 0 | 42.11 |  |  | 0.7731 | no | 0.4047 | BCAS3 |  |
| ERBB3 | auto | ok | 45 | 45 | 18 | 11.11 | 44.44 | 0.7353 | 0.7921 | 0.7731 | no | 0.2417 | ATF1 |  |
| XPO1 | auto | ok | 31 | 31 | 12 | 8.333 | 83.33 | 1 | 1 | 0.7731 | no | 0.25 | USP34 |  |
| ESR1 | auto | ok | 30 | 30 | 13 | 23.08 | 38.46 | 0.8042 | 0.8911 | 0.7731 | no | 0.2834 | LAMA4 |  |
| PAK1 | auto | ok | 25 | 26 | 10 | 10 | 60 | 0.6 | 0.7822 | 0.7731 | no | 0.298 | GDPD4 |  |
| SND1 | auto | ok | 17 | 17 | 17 | 94.12 | 11.76 | 0.8824 | 1 | 0.7731 | no | 0.4286 | BRAF |  |
| CD79A | auto | ok | 13 | 13 | 8 | 0 | 37.5 |  |  | 0.7731 | no | 0.4167 | SHMT1 |  |
| EIF3H | auto | ok | 6 | 6 | 5 | 60 | 40 | 1 | 1 | 0.7731 | no | 0.4286 | RAD21 |  |
| TPR | auto | ok | 5 | 5 | 5 | 40 | 60 | 0.9 | 0.5248 | 0.7731 | no | 0.4386 | NTRK1 |  |
| MAP2K2 | auto | ok | 21 | 21 | 7 | 28.57 | 28.57 | 0.5238 | 0.3267 | 0.7807 | no | 0.2614 | MPND |  |
| ATM | auto | ok | 72 | 72 | 27 | 3.704 | 25.93 | 0.28 | 0.2772 | 0.7992 | no | 0.2396 | C11ORF65 |  |
| SPEN | auto | ok | 61 | 63 | 19 | 0 | 21.05 |  |  | 0.7992 | no | 0.3684 | SPATS1 |  |
| KDM5A | auto | ok | 44 | 45 | 10 | 30 | 50 | 0.119 | 0.2574 | 0.7992 | no | 0.2782 | COLGALT1 |  |
| ABL1 | auto | ok | 21 | 22 | 9 | 0 | 11.11 |  |  | 0.7992 | no | 0.4074 | NUP214 |  |
| FAT1 | auto | ok | 95 | 98 | 11 | 0 | 36.36 |  |  | 0.8196 | no | 0.4545 | SORBS2 |  |
| ELF3 | auto | ok | 47 | 47 | 11 | 27.27 | 45.45 | 0.5952 | 0.1287 | 0.8196 | no | 0.285 | CAPN2 |  |
| CDKN2B | auto | ok | 40 | 40 | 12 | 16.67 | 66.67 | 0.5091 | 0.9307 | 0.8196 | no | 0.3779 | CDKN2A |  |
| BLM | auto | ok | 20 | 20 | 6 | 0 | 33.33 |  |  | 0.8196 | no | 0.4444 | CRTC3 |  |
| MTAP | auto | ok | 8 | 8 | 4 | 25 | 25 | 0.25 | 0.1287 | 0.8196 | no | 0.4544 | CDKN2A |  |
| RAD21 | auto | ok | 38 | 39 | 9 | 44.44 | 55.56 | 0.8333 | 0.495 | 0.8463 | no | 0.2988 | EIF3H |  |
| CTTNBP2 | auto | ok | 5 | 5 | 5 | 0 | 20 |  |  | 0.8513 | no | 0.4922 | MET |  |
| DROSHA | auto | ok | 36 | 36 | 7 | 42.86 | 28.57 | 0.1429 | 0.1485 | 0.8651 | no | 0.3303 | ADAMTS12 |  |
| MGA | auto | ok | 65 | 66 | 15 | 6.667 | 40 | 1 | 1 | 0.9269 | no | 0.2429 | EHD4 |  |
| NOTCH4 | auto | ok | 44 | 46 | 16 | 0 | 75 |  |  | 0.9269 | no | 0.375 | PBX2 |  |
| RICTOR | auto | ok | 34 | 34 | 10 | 30 | 50 | 0.9524 | 0.505 | 0.9269 | no | 0.2855 | FYB1 |  |
| INSR | auto | ok | 31 | 31 | 6 | 16.67 | 50 | 0.5 | 0.8812 | 0.959 | no | 0.2877 | PGPEP1 |  |
| PRKD1 | auto | ok | 16 | 16 | 4 | 50 | 50 | 0.1667 | 0.5941 | 0.959 | no | 0.3295 | FRMD6 |  |
| MPL | auto | ok | 14 | 14 | 6 | 0 | 50 |  |  | 0.959 | no | 0.4444 | LUZP1 |  |
| SMARCB1 | auto | ok | 23 | 23 | 5 | 0 | 20 |  |  | 0.983 | no | 0.4667 | ZNF70 |  |
| GRIN2A | auto | ok | 19 | 19 | 6 | 0 | 16.67 |  |  | 0.983 | no | 0.4444 | ARID1B |  |
| NUF2 | auto | ok | 14 | 14 | 10 | 40 | 10 | 1 | 1 | 0.983 | no | 0.2571 | NOS1AP |  |
| PRKAR1A | auto | ok | 12 | 12 | 6 | 33.33 | 16.67 | 0.4 | 0.6337 | 0.983 | no | 0.2928 | FAM20A |  |
| NF1 | auto | ok | 247 | 253 | 47 | 10.64 | 38.3 | 0.6916 | 0.4653 | 1 | no | 0.2353 | CCDC47 |  |
| ARID1A | auto | ok | 148 | 149 | 31 | 12.9 | 32.26 | 0.8368 | 0.9109 | 1 | no | 0.2296 | SRRM1 |  |
| SMARCA4 | auto | ok | 143 | 146 | 36 | 19.44 | 36.11 | 0.9706 | 0.9703 | 1 | no | 0.2778 | LDLR |  |
| KMT2C | auto | ok | 134 | 136 | 32 | 6.25 | 56.25 | 0.6429 | 0.8614 | 1 | no | 0.2357 | EXOC4 |  |
| KMT2D | auto | ok | 128 | 130 | 47 | 2.128 | 42.55 | 0.5714 | 0.5347 | 1 | no | 0.2357 | ARF3 |  |
| APC | auto | ok | 110 | 112 | 22 | 0 | 45.45 |  |  | 1 | no | 0.4992 | SRP19 |  |
| NOTCH3 | auto | ok | 101 | 102 | 26 | 3.846 | 53.85 | 1 | 1 | 1 | no | 0.248 | BRD4 |  |
| EP300 | auto | ok | 98 | 98 | 22 | 4.545 | 22.73 | 1 | 1 | 1 | no | 0.2338 | MIR1281 |  |
| STK11 | auto | ok | 97 | 98 | 23 | 4.348 | 26.09 | 1 | 1 | 1 | no | 0.2329 | DOCK6 |  |
| KIF5B | auto | ok | 95 | 96 | 96 | 88.54 | 97.92 | 0.2171 | 0.802 | 1 | no | 0.5477 | RET |  |
| ERBB2 | auto | ok | 88 | 92 | 44 | 15.91 | 27.27 | 1 | 1 | 1 | no | 0.224 | RARA |  |
| TSC2 | auto | ok | 88 | 90 | 34 | 17.65 | 23.53 | 0.4563 | 0.8713 | 1 | no | 0.2324 | ZC3H7A |  |
| BRD4 | auto | ok | 86 | 87 | 32 | 3.125 | 50 | 1 | 1 | 1 | no | 0.2411 | AKAP8L |  |
| PTEN | auto | ok | 86 | 86 | 9 | 22.22 | 22.22 | 1 | 1 | 1 | no | 0.2619 | ACTA2 |  |
| DOT1L | auto | ok | 85 | 85 | 22 | 22.73 | 27.27 | 0.8341 | 0.6139 | 1 | no | 0.2605 | PTK2 |  |
| PBRM1 | auto | ok | 85 | 85 | 18 | 5.556 | 61.11 | 1 | 1 | 1 | no | 0.2981 | BAP1 |  |
| ZFHX3 | auto | ok | 85 | 87 | 7 | 0 | 28.57 |  |  | 1 | no | 0.3741 | CFDP1 |  |
| BRCA2 | auto | ok | 79 | 82 | 12 | 8.333 | 41.67 | 1 | 1 | 1 | no | 0.25 | PCDH7 |  |
| POLE | auto | ok | 73 | 74 | 22 | 18.18 | 45.45 | 0.9323 | 0.7624 | 1 | no | 0.2256 | LINC01606 |  |
| ARID2 | auto | ok | 72 | 72 | 13 | 0 | 7.692 |  |  | 1 | no | 0.3846 | CNTN1 |  |
| CTNNB1 | auto | ok | 72 | 72 | 10 | 30 | 20 | 0.5833 | 0.198 | 1 | no | 0.2823 | ULK4 |  |
| KDM6A | auto | ok | 71 | 71 | 8 | 25 | 50 | 0.4 | 0.9703 | 1 | no | 0.2683 | PIGL |  |
| CIC | auto | ok | 68 | 72 | 21 | 0 | 42.86 |  |  | 1 | no | 0.3905 | CSMD1 |  |
| BRCA1 | auto | ok | 66 | 69 | 27 | 7.407 | 55.56 | 0.85 | 0.8515 | 1 | no | 0.2397 | METTL25 |  |
| KMT2A | auto | ok | 65 | 65 | 21 | 19.05 | 23.81 | 0.5417 | 0.7921 | 1 | no | 0.2587 | PHC2 |  |
| FGFR1 | auto | ok | 64 | 64 | 23 | 30.43 | 52.17 | 0.2678 | 0.703 | 1 | no | 0.2384 | IGHMBP2 |  |
| SMAD4 | auto | ok | 62 | 62 | 6 | 0 | 33.33 |  |  | 1 | no | 0.4444 | UNC13C |  |
| PTPRD | auto | ok | 57 | 58 | 8 | 12.5 | 50 | 0.5 | 0.7228 | 1 | no | 0.2729 | LINC01231 |  |
| RTEL1 | auto | ok | 57 | 57 | 24 | 20.83 | 29.17 | 0.4625 | 0.4653 | 1 | no | 0.2619 | CDH4 |  |
| PIK3R1 | auto | ok | 55 | 56 | 3 | 33.33 | 33.33 | 1 | 1 | 1 | no | 0.1818 | LRRFIP1 |  |
| ANKRD11 | auto | ok | 54 | 54 | 6 | 0 | 66.67 |  |  | 1 | no | 0.4444 | FGFR2 |  |
| MDC1 | auto | ok | 52 | 56 | 14 | 7.143 | 71.43 | 0.7143 | 0.6832 | 1 | no | 0.2648 | MUC22 |  |
| FOXA1 | auto | ok | 51 | 52 | 11 | 0 | 18.18 |  |  | 1 | no | 0.5969 | SLC25A21 |  |
| GLI1 | auto | ok | 51 | 53 | 21 | 19.05 | 42.86 | 0.6254 | 0.8218 | 1 | no | 0.2377 | R3HDM2 |  |
| STAT3 | auto | ok | 51 | 51 | 12 | 8.333 | 41.67 | 1 | 1 | 1 | no | 0.25 | CNOT2 |  |
| FANCA | auto | ok | 49 | 49 | 20 | 15 | 40 | 0.807 | 0.802 | 1 | no | 0.2187 | ABHD3 |  |
| ATR | auto | ok | 48 | 48 | 7 | 0 | 57.14 |  |  | 1 | no | 0.4286 | ACP3 |  |
| NF2 | auto | ok | 46 | 48 | 12 | 8.333 | 8.333 | 1 | 1 | 1 | no | 0.3571 | EWSR1 |  |
| SETD2 | auto | ok | 46 | 46 | 17 | 0 | 58.82 |  |  | 1 | no | 0.3725 | TSEN2 |  |
| NSD1 | auto | ok | 45 | 45 | 12 | 8.333 | 33.33 | 0.3636 | 0.7129 | 1 | no | 0.2552 | CHKA |  |
| PTPRS | auto | ok | 44 | 45 | 14 | 7.143 | 50 | 0.5385 | 0.7921 | 1 | no | 0.2704 | PPAN |  |
| PAX8 | auto | ok | 43 | 45 | 10 | 30 | 100 | 1 | 1 | 1 | no | 0.05455 | ACOXL |  |
| PIK3C2G | auto | ok | 43 | 43 | 10 | 30 | 30 | 0.7619 | 0.7822 | 1 | no | 0.261 | GYS2 |  |
| RPTOR | auto | ok | 43 | 43 | 13 | 23.08 | 7.692 | 1 | 1 | 1 | no | 0.2473 | MAP3K3 |  |
| B2M | auto | ok | 42 | 43 | 4 | 0 | 0 |  |  | 1 | no | 0.7386 | TRIM69 |  |
| ERBB4 | auto | ok | 42 | 42 | 5 | 40 | 60 | 0.3 | 0.6535 | 1 | no | 0.3066 | SLC13A3 |  |
| DNMT3B | auto | ok | 40 | 41 | 7 | 0 | 28.57 |  |  | 1 | no | 0.4286 | MSRA |  |
| KEAP1 | auto | ok | 40 | 40 | 8 | 0 | 50 |  |  | 1 | no | 0.4167 | DNAH1 |  |
| NFE2L2 | auto | ok | 40 | 41 | 5 | 0 | 20 |  |  | 1 | no | 0.5653 | ZNF385B |  |
| CARM1 | auto | ok | 39 | 39 | 18 | 22.22 | 22.22 | 0.6729 | 0.9802 | 1 | no | 0.2258 | DNM2 |  |
| EZH2 | auto | ok | 39 | 39 | 14 | 14.29 | 14.29 | 1 | 1 | 1 | no | 0.2755 | CUL1 |  |
| FOXP1 | auto | ok | 39 | 39 | 5 | 20 | 80 | 0.8 | 1 | 1 | no | 0.3 | TMPRSS2 |  |
| NBN | auto | ok | 38 | 38 | 6 | 0 | 66.67 |  |  | 1 | no | 0.4444 | CALB1 |  |
| AXL | auto | ok | 37 | 37 | 10 | 30 | 10 | 1 | 1 | 1 | no | 0.2571 | RNF2 |  |
| PREX2 | auto | ok | 37 | 38 | 4 | 0 | 50 |  |  | 1 | no | 0.5 | TMPRSS2 |  |
| BCOR | auto | ok | 36 | 36 | 7 | 0 | 28.57 |  |  | 1 | no | 0.4286 | KMT2D |  |
| IGF1R | auto | ok | 36 | 37 | 12 | 8.333 | 16.67 | 1 | 1 | 1 | no | 0.25 | ACBD6 |  |
| TBX3 | auto | ok | 36 | 36 | 6 | 0 | 50 |  |  | 1 | no | 0.4444 | ARID1A |  |
| ASXL1 | auto | ok | 35 | 35 | 9 | 0 | 66.67 |  |  | 1 | no | 0.4181 | DNMT3B |  |
| KDM5C | auto | ok | 33 | 34 | 10 | 30 | 50 | 0.9524 | 0.4653 | 1 | no | 0.269 | KANTR |  |
| PTPRT | auto | ok | 33 | 33 | 6 | 33.33 | 50 | 0.9 | 0.9703 | 1 | no | 0.2862 | ZHX3 |  |
| RUNX1 | auto | ok | 33 | 34 | 6 | 16.67 | 33.33 | 0.3333 | 0.9406 | 1 | no | 0.2867 | IRF2BP2 |  |
| BCL2L11 | auto | ok | 32 | 32 | 14 | 21.43 | 14.29 | 1 | 1 | 1 | no | 0.2936 | ACOXL |  |
| DNMT3A | auto | ok | 32 | 32 | 11 | 36.36 | 45.45 | 0.9545 | 1 | 1 | no | 0.345 | DTNB |  |
| JAK3 | auto | ok | 32 | 34 | 9 | 0 | 55.56 |  |  | 1 | no | 0.4074 | CIMAP1D |  |
| PIK3R2 | auto | ok | 32 | 32 | 10 | 10 | 20 | 1 | 1 | 1 | no | 0.2571 | PRKACA |  |
| RPS6KB2 | auto | ok | 32 | 32 | 9 | 33.33 | 11.11 | 1 | 1 | 1 | no | 0.2619 | ANO5 |  |
| ARID5B | auto | ok | 31 | 33 | 6 | 0 | 50 |  |  | 1 | no | 0.4444 | KIF5B |  |
| CSDE1 | auto | ok | 31 | 31 | 5 | 0 | 60 |  |  | 1 | no | 0.4667 | LRRC8B |  |
| ERCC2 | auto | ok | 31 | 32 | 8 | 25 | 87.5 | 0.75 | 0.8614 | 1 | no | 0.2702 | PPP1R13L |  |
| PDGFRA | auto | ok | 31 | 33 | 9 | 55.56 | 44.44 | 0.9921 | 0.3366 | 1 | no | 0.2916 | EXOC1 |  |
| TP53BP1 | auto | ok | 31 | 31 | 11 | 18.18 | 54.55 | 0.9167 | 0.8911 | 1 | no | 0.3282 | PPIP5K1 |  |
| EZH1 | auto | ok | 30 | 31 | 11 | 0 | 18.18 |  |  | 1 | no | 0.3939 | RAB37 |  |
| GNAS | auto | ok | 30 | 30 | 8 | 25 | 0 | 1 | 1 | 1 | no | 0.2679 | PTGIS |  |
| AXIN1 | auto | ok | 29 | 29 | 6 | 50 | 83.33 | 0.5 | 0.6733 | 1 | no | 0.3402 | NPRL3 |  |
| CALR | auto | ok | 29 | 29 | 9 | 44.44 | 0 | 1 | 1 | 1 | no | 0.2619 | FARSA |  |
| CBL | auto | ok | 29 | 29 | 4 | 25 | 25 | 1 | 1 | 1 | no | 0.3817 | KMT2A |  |
| FUBP1 | auto | ok | 29 | 29 | 11 | 18.18 | 45.45 | 0.8333 | 0.9802 | 1 | no | 0.2536 | ATP12A |  |
| JAK1 | auto | ok | 29 | 30 | 8 | 12.5 | 12.5 | 1 | 1 | 1 | no | 0.3214 | GLIS1 |  |
| NTRK2 | auto | ok | 29 | 30 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | COL5A1 |  |
| AXIN2 | auto | ok | 28 | 28 | 4 | 0 | 75 |  |  | 1 | no | 0.5 | PLXDC1 |  |
| KDR | auto | ok | 28 | 29 | 9 | 11.11 | 55.56 | 0.5556 | 0.2673 | 1 | no | 0.2824 | PAICS |  |
| RARA | auto | ok | 28 | 28 | 14 | 28.57 | 28.57 | 1 | 1 | 1 | no | 0.2449 | KRT36 |  |
| TEK | auto | ok | 28 | 28 | 9 | 0 | 22.22 |  |  | 1 | no | 0.4074 | TSPAN9 |  |
| TET1 | auto | ok | 28 | 28 | 7 | 14.29 | 14.29 | 1 | 1 | 1 | no | 0.2755 | CAMTA1 |  |
| IKBKE | auto | ok | 27 | 28 | 12 | 25 | 8.333 | 0.2727 | 0.6139 | 1 | no | 0.2576 | RASAL2 |  |
| MAP2K4 | auto | ok | 27 | 27 | 7 | 14.29 | 0 | 1 | 1 | 1 | no | 0.398 | DNAH9 |  |
| PPM1D | auto | ok | 27 | 29 | 12 | 8.333 | 25 | 1 | 1 | 1 | no | 0.2502 | BCAS3 |  |
| SMAD3 | auto | ok | 27 | 27 | 9 | 11.11 | 55.56 | 1 | 1 | 1 | no | 0.3051 | IQCH |  |
| TRAF7 | auto | ok | 27 | 27 | 10 | 10 | 10 | 1 | 1 | 1 | no | 0.2238 | E4F1 |  |
| IRS2 | auto | ok | 26 | 26 | 9 | 11.11 | 77.78 | 0.7778 | 0.9703 | 1 | no | 0.2624 | DOCK9 |  |
| SOX9 | auto | ok | 26 | 27 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | GRB2 |  |
| UPF1 | auto | ok | 26 | 26 | 6 | 33.33 | 33.33 | 0.7 | 0.2079 | 1 | no | 0.3101 | CEP89 |  |
| CARD11 | auto | ok | 25 | 25 | 7 | 0 | 57.14 |  |  | 1 | no | 0.4286 | CREB5 |  |
| FLT1 | auto | ok | 25 | 25 | 3 | 33.33 | 33.33 | 0.3333 | 0.495 | 1 | no | 0.1957 | IRS2 |  |
| MALT1 | auto | ok | 25 | 25 | 7 | 42.86 | 14.29 | 1 | 1 | 1 | no | 0.2755 | SEC11C |  |
| ASXL2 | auto | ok | 24 | 24 | 7 | 28.57 | 57.14 | 0.2857 | 0.495 | 1 | no | 0.3476 | PPM1G |  |
| DIS3 | auto | ok | 24 | 24 | 6 | 16.67 | 16.67 | 1 | 1 | 1 | no | 0.4286 | NALF1 |  |
| TP63 | auto | ok | 24 | 24 | 7 | 28.57 | 14.29 | 1 | 1 | 1 | no | 0.2755 | LINC00578 |  |
| AKT1 | auto | ok | 23 | 23 | 7 | 0 | 57.14 |  |  | 1 | no | 0.4523 | LAMTOR1 |  |
| MAP3K13 | auto | ok | 23 | 24 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | MASP1 |  |
| MDM2 | auto | ok | 23 | 23 | 17 | 23.53 | 41.18 | 0.9118 | 0.8812 | 1 | no | 0.2667 | CADM3 |  |
| MED12 | auto | ok | 23 | 23 | 6 | 0 | 33.33 |  |  | 1 | no | 0.4444 | LINC01281 |  |
| RAF1 | auto | ok | 23 | 23 | 12 | 25 | 58.33 | 0.9545 | 1 | 1 | no | 0.25 | PPP4R2 |  |
| TCF7L2 | auto | ok | 23 | 23 | 4 | 0 | 25 |  |  | 1 | no | 0.5 | FAM204A |  |
| EIF4A2 | auto | ok | 22 | 22 | 6 | 0 | 16.67 |  |  | 1 | no | 0.5449 | RFC4 |  |
| PARP1 | auto | ok | 22 | 22 | 10 | 10 | 70 | 0.7 | 0.9307 | 1 | no | 0.3189 | LIN9 |  |
| PIK3CD | auto | ok | 22 | 23 | 8 | 12.5 | 50 | 0.5 | 0.9109 | 1 | no | 0.2693 | MARK4 |  |
| PLK2 | auto | ok | 22 | 22 | 7 | 28.57 | 57.14 | 0.2857 | 0.3366 | 1 | no | 0.3391 | PDE4D |  |
| PPP2R1A | auto | ok | 22 | 23 | 6 | 50 | 33.33 | 0.8 | 0.7525 | 1 | no | 0.2901 | FRMD4A |  |
| PRKCI | auto | ok | 22 | 22 | 6 | 33.33 | 33.33 | 0.6 | 0.8218 | 1 | no | 0.2145 | VPS8 |  |
| RAD54L | auto | ok | 22 | 22 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | XPO7 |  |
| FLT3 | auto | ok | 21 | 21 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | unknown |  |
| RPS6KA4 | auto | ok | 21 | 22 | 13 | 38.46 | 53.85 | 0.2475 | 0.5743 | 1 | no | 0.2559 | MACROD1 |  |
| SMARCD1 | auto | ok | 21 | 21 | 8 | 12.5 | 50 | 0.5714 | 0.7921 | 1 | no | 0.2771 | SPATS2 |  |
| ASPSCR1 | auto | ok | 20 | 20 | 20 | 70 | 20 | 0.9391 | 0.8317 | 1 | no | 0.5488 | TFE3 |  |
| FLCN | auto | ok | 20 | 21 | 10 | 50 | 60 | 0.881 | 0.9505 | 1 | no | 0.2722 | COPS3 |  |
| MSH6 | auto | ok | 20 | 20 | 6 | 0 | 50 |  |  | 1 | no | 0.4444 | RHOQ |  |
| PIK3CB | auto | ok | 20 | 20 | 5 | 0 | 60 |  |  | 1 | no | 0.4667 | FOXL2 |  |
| PIK3R3 | auto | ok | 20 | 20 | 3 | 33.33 | 33.33 | 1 | 1 | 1 | no | 0.1818 | ARHGEF10L |  |
| SMYD3 | auto | ok | 20 | 20 | 4 | 25 | 50 | 0.6667 | 0.2178 | 1 | no | 0.1664 | COQ5 |  |
| SRC | auto | ok | 20 | 20 | 4 | 25 | 0 | 1 | 1 | 1 | no | 0.3214 | RPN2 |  |
| CASP8 | auto | ok | 19 | 19 | 9 | 22.22 | 33.33 | 1 | 1 | 1 | no | 0.2619 | ORC2 |  |
| CHEK2 | auto | ok | 19 | 19 | 6 | 16.67 | 66.67 | 0.6667 | 0.3564 | 1 | no | 0.4862 | TTC28 |  |
| CSF1R | auto | ok | 19 | 19 | 9 | 0 | 55.56 |  |  | 1 | no | 0.4074 | MAN2A1 |  |
| LATS1 | auto | ok | 19 | 19 | 6 | 16.67 | 66.67 | 0.6667 | 0.5149 | 1 | no | 0.296 | MACROD2 |  |
| PTCH1 | auto | ok | 19 | 19 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | STRN |  |
| YAP1 | auto | ok | 19 | 19 | 5 | 20 | 20 | 1 | 1 | 1 | no | 0.3 | MAML2 |  |
| FGFR4 | auto | ok | 18 | 18 | 9 | 11.11 | 44.44 | 1 | 1 | 1 | no | 0.2619 | ZNF346 |  |
| PGR | auto | ok | 18 | 18 | 7 | 14.29 | 71.43 | 0.7143 | 0.9802 | 1 | no | 0.2758 | CCND1 |  |
| VTCN1 | auto | ok | 18 | 18 | 5 | 20 | 60 | 0.6 | 0.9901 | 1 | no | 0.3002 | VSIG10 |  |
| CSF3R | auto | ok | 17 | 18 | 6 | 33.33 | 66.67 | 0.4 | 0.5941 | 1 | no | 0.2938 | TFAP2E |  |
| EGFL7 | auto | ok | 17 | 17 | 4 | 50 | 25 | 1 | 1 | 1 | no | 0.4286 | PNPLA7 |  |
| ERCC5 | auto | ok | 17 | 17 | 8 | 0 | 37.5 |  |  | 1 | no | 0.5 | USH2A |  |
| LDLR | auto | ok | 17 | 17 | 15 | 20 | 60 | 0.7253 | 0.8515 | 1 | no | 0.4455 | SMARCA4 |  |
| PLCG2 | auto | ok | 17 | 17 | 6 | 33.33 | 33.33 | 1 | 1 | 1 | no | 0.2857 | WWOX |  |
| SDC4 | auto | ok | 17 | 17 | 17 | 82.35 | 0 | 1 | 1 | 1 | no | 0.4286 | ROS1 |  |
| SEC16A | auto | ok | 17 | 17 | 17 | 23.53 | 47.06 | 1 | 1 | 1 | no | 0.4286 | NOTCH1 |  |
| TAP1 | auto | ok | 17 | 17 | 7 | 0 | 57.14 |  |  | 1 | no | 0.4286 | HLA-DMB |  |
| TGFBR1 | auto | ok | 17 | 17 | 4 | 25 | 75 | 0.75 | 0.7525 | 1 | no | 0.3583 | COL15A1 |  |
| CTCF | auto | ok | 16 | 16 | 4 | 0 | 25 |  |  | 1 | no | 0.5 | CARMIL2 |  |
| CUL1 | auto | ok | 16 | 16 | 14 | 35.71 | 0 | 1 | 1 | 1 | no | 0.3321 | BRAF |  |
| DICER1 | auto | ok | 16 | 16 | 5 | 0 | 40 |  |  | 1 | no | 0.5813 | DLGAP5 |  |
| PDGFRB | auto | ok | 16 | 16 | 5 | 40 | 60 | 0.9 | 0.9505 | 1 | no | 0.3008 | SH3PXD2B |  |
| RAD51B | auto | ok | 16 | 16 | 2 | 50 | 100 | 1 | 1 | 1 | no | 0.2727 | CCND2 |  |
| RECQL | auto | ok | 16 | 16 | 4 | 0 | 50 |  |  | 1 | no | 0.5 | TMEM117 |  |
| SMO | auto | ok | 16 | 16 | 8 | 0 | 0 |  |  | 1 | no | 0.4167 | MKRN1 |  |
| SUZ12 | auto | ok | 16 | 16 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | KMT2A |  |
| CCND1 | auto | ok | 15 | 15 | 7 | 0 | 57.14 |  |  | 1 | no | 0.4286 | PICALM |  |
| CCNE1 | auto | ok | 15 | 15 | 5 | 20 | 20 | 1 | 1 | 1 | no | 0.3 | CLEC1A |  |
| TNFRSF14 | auto | ok | 15 | 15 | 4 | 25 | 0 | 1 | 1 | 1 | no | 0.2727 | PLCH2 |  |
| CDK4 | auto | ok | 14 | 14 | 6 | 33.33 | 16.67 | 1 | 1 | 1 | no | 0.2857 | LGR5 |  |
| E2F3 | auto | ok | 14 | 14 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | FGD2 |  |
| NUP93 | auto | ok | 14 | 14 | 5 | 0 | 0 |  |  | 1 | no | 0.4667 | MT1H |  |
| RIT1 | auto | ok | 14 | 14 | 6 | 0 | 16.67 |  |  | 1 | no | 0.4444 | KCNN3 |  |
| SESN3 | auto | ok | 14 | 15 | 5 | 20 | 0 | 1 | 1 | 1 | no | 0.1091 | GRIA4 |  |
| SF3B1 | auto | ok | 14 | 15 | 5 | 20 | 40 | 1 | 1 | 1 | no | 0.3 | PLEKHM3 |  |
| MST1R | auto | ok | 13 | 13 | 5 | 20 | 60 | 0.6 | 0.6733 | 1 | no | 0.3061 | IP6K1 |  |
| NPM1 | auto | ok | 13 | 13 | 4 | 50 | 25 | 1 | 1 | 1 | no | 0.3214 | FGF18 |  |
| YES1 | auto | ok | 13 | 13 | 6 | 16.67 | 16.67 | 1 | 1 | 1 | no | 0.4286 | GRK3 |  |
| AKT3 | auto | ok | 12 | 12 | 3 | 33.33 | 66.67 | 1 | 1 | 1 | no | 0.3636 | SDCCAG8 |  |
| CEBPA | auto | ok | 12 | 12 | 6 | 0 | 16.67 |  |  | 1 | no | 0.4444 | UBA2 |  |
| DUSP4 | auto | ok | 12 | 12 | 6 | 0 | 33.33 |  |  | 1 | no | 0.4444 | KIF13B |  |
| FYN | auto | ok | 12 | 12 | 3 | 66.67 | 66.67 | 1 | 0.9802 | 1 | no | 0.1822 | LAMA2 |  |
| IDH1 | auto | ok | 12 | 12 | 4 | 25 | 0 | 1 | 1 | 1 | no | 0.2727 | CMKLR2-AS |  |
| PAX5 | auto | ok | 12 | 12 | 3 | 33.33 | 66.67 | 0.6667 | 0.9802 | 1 | no | 0.1822 | ADAMTSL1 |  |
| PIK3C3 | auto | ok | 12 | 12 | 1 | 100 | 0 | 1 | 1 | 1 | no | 0.5455 | SLC35D4 |  |
| PIK3CA | auto | ok | 12 | 12 | 5 | 0 | 20 |  |  | 1 | no | 0.5945 | NAALADL2 |  |
| SDHA | auto | ok | 12 | 12 | 5 | 40 | 40 | 0.7 | 0.4851 | 1 | no | 0.3112 | AHRR |  |
| ERCC4 | auto | ok | 11 | 11 | 5 | 20 | 20 | 1 | 1 | 1 | no | 0.3 | RBFOX1 |  |
| FGF19 | auto | ok | 11 | 11 | 8 | 0 | 75 |  |  | 1 | no | 0.4167 | CCND1 |  |
| GSK3B | auto | ok | 11 | 11 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | FSTL1 |  |
| MKRN1 | auto | ok | 11 | 11 | 11 | 63.64 | 0 | 1 | 1 | 1 | no | 0.4959 | BRAF |  |
| MLH1 | auto | ok | 11 | 14 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | SLC4A7 |  |
| PHF7 | auto | ok | 11 | 11 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | BAP1 |  |
| RAD51D | auto | ok | 11 | 11 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | ASIC2 |  |
| RELA | auto | ok | 11 | 11 | 6 | 33.33 | 100 | 1 | 1 | 1 | no | 0.4545 | ZFTA |  |
| SOS1 | auto | ok | 11 | 11 | 3 | 33.33 | 0 | 1 | 1 | 1 | no | 0.1818 | ALK |  |
| TACC2 | auto | ok | 11 | 11 | 8 | 62.5 | 100 | 1 | 1 | 1 | no | 0.5455 | FGFR2 |  |
| BMPR1A | auto | ok | 10 | 10 | 1 | 100 | 0 | 1 | 1 | 1 | no | 0.5455 | DACH2 |  |
| CDK8 | auto | ok | 10 | 10 | 4 | 25 | 0 | 1 | 1 | 1 | no | 0.1364 | CACNA1I |  |
| CUL3 | auto | ok | 10 | 10 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | CERKL |  |
| EPAS1 | auto | ok | 10 | 10 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | SBF2 |  |
| IKZF3 | auto | ok | 10 | 10 | 8 | 25 | 100 | 1 | 1 | 1 | no | 0.3409 | ERBB2 |  |
| MDM4 | auto | ok | 10 | 10 | 4 | 25 | 25 | 1 | 1 | 1 | no | 0.3214 | EFCAB6 |  |
| RXRA | auto | ok | 10 | 10 | 3 | 33.33 | 66.67 | 1 | 1 | 1 | no | 0.1818 | RAPGEF1 |  |
| TRIM24 | auto | ok | 10 | 10 | 9 | 88.89 | 0 | 1 | 1 | 1 | no | 0.4848 | BRAF |  |
| ARAF | auto | ok | 9 | 10 | 2 | 50 | 50 | 1 | 1 | 1 | no | 0.2727 | SCML2 |  |
| LRP1 | auto | ok | 9 | 9 | 5 | 0 | 20 |  |  | 1 | no | 0.7038 | NAB2 |  |
| NTHL1 | auto | ok | 9 | 9 | 3 | 33.33 | 66.67 | 0.6667 | 0.5347 | 1 | no | 0.1942 | ABCA3 |  |
| PPARG | auto | ok | 9 | 9 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | IQSEC1 |  |
| RAD51 | auto | ok | 9 | 9 | 3 | 66.67 | 0 | 1 | 1 | 1 | no | 0.3636 | EGFR |  |
| SESN1 | auto | ok | 9 | 9 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | PPIL6 |  |
| STK40 | auto | ok | 9 | 9 | 3 | 33.33 | 0 | 1 | 1 | 1 | no | 0.1818 | ASXL2 |  |
| ACOXL | auto | ok | 8 | 8 | 8 | 12.5 | 50 | 0.5 | 0.505 | 1 | no | 0.3856 | BCL2L11 |  |
| ETV4 | auto | ok | 8 | 8 | 6 | 33.33 | 16.67 | 1 | 1 | 1 | no | 0.396 | TMPRSS2 |  |
| HNF1A | auto | ok | 8 | 8 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | CFAP251 |  |
| IKZF1 | auto | ok | 8 | 8 | 1 | 100 | 0 | 1 | 1 | 1 | no | 0.5455 | BBS9 |  |
| IRF4 | auto | ok | 8 | 8 | 5 | 0 | 20 |  |  | 1 | no | 0.4667 | IFNGR1 |  |
| KRAS | auto | ok | 8 | 8 | 3 | 33.33 | 0 | 1 | 1 | 1 | no | 0.1818 | MS4A2 |  |
| MAPK1 | auto | ok | 8 | 8 | 2 | 50 | 50 | 0.5 | 0.5446 | 1 | no | 0.2847 | GNL2 |  |
| MCL1 | auto | ok | 8 | 8 | 6 | 0 | 50 |  |  | 1 | no | 0.4444 | NOS1AP |  |
| MSI1 | auto | ok | 8 | 8 | 5 | 20 | 40 | 1 | 1 | 1 | no | 0.3108 | GCN1 |  |
| MYD88 | auto | ok | 8 | 8 | 6 | 33.33 | 33.33 | 0.6 | 0.7129 | 1 | no | 0.291 | ACAA1 |  |
| PDE4D | auto | ok | 8 | 8 | 6 | 16.67 | 66.67 | 0.6667 | 0.2277 | 1 | no | 0.3801 | MAP3K1 |  |
| PTK2 | auto | ok | 8 | 8 | 8 | 25 | 50 | 1 | 1 | 1 | no | 0.3862 | AGO2 |  |
| SLC25A21 | auto | ok | 8 | 8 | 8 | 0 | 12.5 |  |  | 1 | no | 0.5 | FOXA1 |  |
| CACNA1A | auto | ok | 7 | 7 | 7 | 14.29 | 42.86 | 0.5 | 0.5248 | 1 | no | 0.2855 | DNMT1 |  |
| CDK6 | auto | ok | 7 | 7 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | PPP1R9A |  |
| CRKL | auto | ok | 7 | 7 | 4 | 0 | 50 |  |  | 1 | no | 0.5 | RASGEF1C |  |
| DCUN1D1 | auto | ok | 7 | 7 | 3 | 66.67 | 0 | 1 | 1 | 1 | no | 0.1818 | ASCC3 |  |
| DGKH | auto | ok | 7 | 7 | 7 | 0 | 28.57 |  |  | 1 | no | 0.4762 | FH |  |
| DNM2 | auto | ok | 7 | 7 | 7 | 28.57 | 42.86 | 0.7143 | 0.9406 | 1 | no | 0.2298 | CARM1 |  |
| GAB2 | auto | ok | 7 | 7 | 5 | 20 | 20 | 1 | 1 | 1 | no | 0.3 | FGFR2 |  |
| IGF2 | auto | ok | 7 | 8 | 4 | 50 | 75 | 0.5 | 1 | 1 | no | 0.4286 | PHF21A |  |
| KREMEN1 | auto | ok | 7 | 7 | 7 | 28.57 | 42.86 | 1 | 1 | 1 | no | 0.2755 | NF2 |  |
| NR4A3 | auto | ok | 7 | 7 | 7 | 57.14 | 71.43 | 0.7143 | 1 | 1 | no | 0.4286 | EWSR1 |  |
| ST7 | auto | ok | 7 | 7 | 7 | 14.29 | 0 | 1 | 1 | 1 | no | 0.5725 | MET |  |
| TACC1 | auto | ok | 7 | 7 | 4 | 25 | 100 | 1 | 1 | 1 | no | 0.2727 | FGFR1 |  |
| USP34 | auto | ok | 7 | 7 | 7 | 14.29 | 85.71 | 1 | 1 | 1 | no | 0.4286 | XPO1 |  |
| CAMTA1 | auto | ok | 6 | 6 | 5 | 20 | 20 | 1 | 1 | 1 | no | 0.3 | TET1 |  |
| CDKN2C | auto | ok | 6 | 6 | 3 | 33.33 | 66.67 | 0.6667 | 0.3267 | 1 | no | 0.2039 | NKAIN1 |  |
| CPM | auto | ok | 6 | 6 | 5 | 40 | 20 | 1 | 1 | 1 | no | 0.3 | CSF1R |  |
| CTNNA3 | auto | ok | 6 | 6 | 5 | 60 | 20 | 0.6 | 0.2277 | 1 | no | 0.3973 | FGFR2 |  |
| DTNB | auto | ok | 6 | 6 | 6 | 16.67 | 0 | 1 | 1 | 1 | no | 0.3571 | DNMT3A |  |
| EPHX3 | auto | ok | 6 | 6 | 6 | 16.67 | 33.33 | 1 | 1 | 1 | no | 0.2317 | BRD4 |  |
| GALNT11 | auto | ok | 6 | 6 | 6 | 0 | 66.67 |  |  | 1 | no | 0.6667 | KMT2C |  |
| GRB7 | auto | ok | 6 | 6 | 5 | 0 | 20 |  |  | 1 | no | 0.5333 | ERBB2 |  |
| NOS1AP | auto | ok | 6 | 6 | 6 | 0 | 0 |  |  | 1 | no | 0.4444 | FLCN |  |
| SMARCA2 | auto | ok | 6 | 6 | 3 | 33.33 | 66.67 | 0.6667 | 0.8218 | 1 | no | 0.1857 | EWSR1 |  |
| SOX2 | auto | ok | 6 | 6 | 4 | 0 | 50 |  |  | 1 | no | 0.4483 | SOX2-OT |  |
| TANC2 | auto | ok | 6 | 6 | 6 | 16.67 | 33.33 | 1 | 1 | 1 | no | 0.2857 | RNF43 |  |
| TFG | auto | ok | 6 | 6 | 6 | 100 | 100 | 1 | 1 | 1 | no | 0.2727 | ROS1 |  |
| ATE1 | auto | ok | 5 | 5 | 3 | 33.33 | 0 | 1 | 1 | 1 | no | 0.5455 | FGFR2 |  |
| BCL2L1 | auto | ok | 5 | 5 | 2 | 50 | 50 | 0.5 | 0.495 | 1 | no | 0.2866 | FSIP1 |  |
| CNTNAP2 | auto | ok | 5 | 5 | 3 | 66.67 | 0 | 1 | 1 | 1 | no | 0.3636 | MET |  |
| CREM | auto | ok | 5 | 5 | 5 | 20 | 100 | 1 | 1 | 1 | no | 0.3273 | EWSR1 |  |
| EXOC4 | auto | ok | 5 | 5 | 5 | 0 | 20 |  |  | 1 | no | 0.6 | KMT2C |  |
| FBXL7 | auto | ok | 5 | 5 | 4 | 50 | 75 | 0.5 | 0.5941 | 1 | no | 0.4366 | RET |  |
| GPATCH8 | auto | ok | 5 | 5 | 4 | 25 | 50 | 0.5 | 0.8515 | 1 | no | 0.3239 | MAP3K14 |  |
| HCN1 | auto | ok | 5 | 5 | 1 | 100 | 0 | 1 | 1 | 1 | no | 0.5455 | DROSHA |  |
| LRP1B | auto | ok | 5 | 5 | 2 | 50 | 50 | 1 | 1 | 1 | no | 0.2727 | ALK |  |
| MAD1L1 | auto | ok | 5 | 5 | 5 | 40 | 0 | 1 | 1 | 1 | no | 0.3857 | TMPRSS2 |  |
| PICALM | auto | ok | 5 | 5 | 5 | 0 | 80 |  |  | 1 | no | 0.4667 | CCND1 |  |
| PRKAG2 | auto | ok | 5 | 5 | 5 | 20 | 40 | 1 | 1 | 1 | no | 0.3429 | KMT2C |  |
| RAC2 | auto | ok | 5 | 5 | 2 | 50 | 0 | 1 | 1 | 1 | no | 0.2727 | ATP6V0A4 |  |
| SBNO2 | auto | ok | 5 | 5 | 4 | 0 | 75 |  |  | 1 | no | 0.5 | TCF3 |  |
| SDHAF2 | auto | ok | 5 | 6 | 1 | 100 | 0 | 1 | 1 | 1 | no | 0.5455 | MGMT |  |
| SFPQ | auto | ok | 5 | 5 | 4 | 25 | 50 | 1 | 1 | 1 | no | 0.4286 | TFE3 |  |
| SHANK2 | auto | ok | 5 | 5 | 4 | 0 | 50 |  |  | 1 | no | 0.5 | FGF3 |  |
| SMARCE1 | auto | ok | 5 | 5 | 1 | 100 | 100 | 1 | 1 | 1 | no | 0.5455 | DNAJC17 |  |
| SMURF2 | auto | ok | 5 | 5 | 4 | 25 | 75 | 0.75 | 0.3564 | 1 | no | 0.3374 | CD79B |  |
| SPTBN1 | auto | ok | 5 | 5 | 4 | 50 | 0 | 1 | 1 | 1 | no | 0.5455 | ALK |  |
| THSD4 | auto | ok | 5 | 5 | 4 | 25 | 25 | 1 | 1 | 1 | no | 0.3214 | STAT5A |  |
| TMEM117 | auto | ok | 5 | 5 | 5 | 0 | 0 |  |  | 1 | no | 0.4667 | RECQL |  |
| WWOX | auto | ok | 5 | 5 | 5 | 0 | 40 |  |  | 1 | no | 0.4667 | TP53 |  |
| CDKN2A | auto | ok | 185 | 187 | 21 | 4.762 | 0 |  |  |  | no | 0.2857 | CDKN2B |  |
| KMT2B | unresolved | failed | 68 | 69 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| NSD3 | auto | ok | 58 | 58 | 28 | 14.29 | 0 |  |  |  | no | 0.1429 | FGFR1 |  |
| RECQL4 | unresolved | failed | 54 | 54 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| PRKN | auto | ok | 40 | 40 | 4 | 50 | 0 |  |  |  | no | 0.25 | ESR1 |  |
| NSD2 | auto | ok | 37 | 37 | 9 | 33.33 | 0 |  |  |  | no | 0.3333 | FGFR3 |  |
| MAP3K14 | unresolved | failed | 36 | 36 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| TCF3 | unresolved | failed | 35 | 36 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| RASA1 | auto | ok | 32 | 32 | 4 | 0 | 50 |  |  |  | no | 0.5 | KIAA0825 |  |
| STAG2 | auto | ok | 31 | 31 | 2 | 0 | 0 |  |  |  | no | 0.5 | KCNJ15 |  |
| AR | auto | ok | 30 | 30 | 2 | 0 | 50 |  |  |  | no | 0.5 | CAGE1 |  |
| CDC73 | auto | ok | 28 | 28 | 2 | 0 | 50 |  |  |  | no | 0.5 | HIVEP3 |  |
| JAK2 | auto | ok | 27 | 27 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | CDC37L1 |  |
| CDKN1B | auto | ok | 24 | 24 | 3 | 0 | 0 |  |  |  | no | 0.3333 | BCL2L14 |  |
| MSH2 | auto | ok | 23 | 23 | 4 | 0 | 0 |  |  |  | no | 0.25 | CTNNA2 |  |
| NFKBIA | auto | ok | 23 | 23 | 1 | 0 | 100 |  |  |  | no | 1 | CHCT1 |  |
| MEN1 | auto | ok | 22 | 22 | 3 | 0 | 0 |  |  |  | no | 0.3333 | MAML2 |  |
| MYC | auto | ok | 22 | 22 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | BCR |  |
| TET2 | auto | ok | 21 | 23 | 3 | 0 | 66.67 |  |  |  | no | 0.6667 | ARHGEF38 |  |
| ERF | unresolved | failed | 20 | 20 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| KIT | auto | ok | 20 | 20 | 2 | 0 | 0 |  |  |  | no | 0.5 | KDR |  |
| CDKN2B-AS1 | unresolved | failed | 19 | 19 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| EPHA3 | auto | ok | 19 | 19 | 2 | 0 | 0 |  |  |  | no | 0.5 | MAPK6 |  |
| NEGR1 | auto | ok | 19 | 19 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | HELZ |  |
| PALB2 | auto | ok | 19 | 19 | 5 | 0 | 0 |  |  |  | no | 0.2 | ERN2 |  |
| SEPTIN14 | unresolved | failed | 19 | 19 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| BBC3 | auto | ok | 18 | 18 | 8 | 50 | 0 |  |  |  | no | 0.375 | SAE1 |  |
| GATA3 | auto | ok | 18 | 18 | 2 | 0 | 0 |  |  |  | no | 0.5 | CAMK1D |  |
| AMER1 | auto | ok | 17 | 17 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | ASB12 |  |
| COP1 | auto | ok | 17 | 17 | 2 | 50 | 0 |  |  |  | no | 0.5 | DNAJB4 |  |
| DDR2 | auto | ok | 17 | 17 | 5 | 0 | 0 |  |  |  | no | 0.2 | CD247 |  |
| EPHA5 | auto | ok | 17 | 17 | 3 | 0 | 33.33 |  |  |  | no | 0.6667 | TECRL |  |
| INPP4B | auto | ok | 17 | 18 | 1 | 100 | 0 |  |  |  | no | 1 | FREM3 |  |
| RAD50 | auto | ok | 17 | 18 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | FBN2 |  |
| BCL6 | auto | ok | 16 | 16 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | DGKG |  |
| ERRFI1 | auto | ok | 16 | 16 | 3 | 0 | 0 |  |  |  | no | 0.6667 | ZNF135 |  |
| LATS2 | auto | ok | 16 | 16 | 1 | 0 | 0 |  |  |  | no | 1 | ATP8A2 |  |
| MITF | auto | ok | 16 | 16 | 1 | 0 | 100 |  |  |  | no | 1 | FOXP1 |  |
| MSI2 | auto | ok | 16 | 16 | 2 | 0 | 50 |  |  |  | no | 0.5 | BCAS3 |  |
| NCOA4 | auto | ok | 16 | 16 | 16 | 25 | 0 |  |  |  | no | 1 | RET |  |
| RAD52 | auto | ok | 16 | 16 | 1 | 0 | 100 |  |  |  | no | 1 | IRAK2 |  |
| TGFBR2 | auto | ok | 16 | 16 | 2 | 0 | 50 |  |  |  | no | 0.5 | GADL1 |  |
| TRAP1 | unresolved | failed | 16 | 16 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| CREB3L1 | unresolved | failed | 15 | 15 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| INPP4A | auto | ok | 15 | 15 | 7 | 14.29 | 0 |  |  |  | no | 0.2857 | MGAT4A |  |
| MSH3 | auto | ok | 15 | 16 | 4 | 0 | 0 |  |  |  | no | 0.5 | TERT |  |
| GPS2 | auto | ok | 14 | 16 | 2 | 0 | 0 |  |  |  | no | 0.5 | DNAH9 |  |
| HGF | auto | ok | 14 | 14 | 4 | 0 | 50 |  |  |  | no | 0.5 | AATF |  |
| LINC00114 | unresolved | failed | 14 | 14 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| TAP2 | auto | ok | 14 | 14 | 2 | 0 | 50 |  |  |  | no | 0.5 | ACLY |  |
| TRAF2 | unresolved | failed | 14 | 14 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| CYLD | auto | ok | 13 | 13 | 1 | 0 | 0 |  |  |  | no | 1 | TACC1 |  |
| DAXX | auto | ok | 13 | 13 | 1 | 0 | 100 |  |  |  | no | 1 | TFAP2B |  |
| MAP2K1 | auto | ok | 13 | 13 | 4 | 0 | 25 |  |  |  | no | 0.25 | OS9 |  |
| MEF2B | unresolved | failed | 13 | 13 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| MRE11 | auto | ok | 13 | 13 | 5 | 0 | 20 |  |  |  | no | 0.2 | ATM |  |
| PAK5 | auto | ok | 13 | 13 | 2 | 50 | 0 |  |  |  | no | 0.5 | KDM5C |  |
| RAD51C | auto | ok | 13 | 13 | 6 | 0 | 0 |  |  |  | no | 0.1667 | BCAS3 |  |
| SPOP | auto | ok | 13 | 13 | 2 | 0 | 100 |  |  |  | no | 0.5 | GDPD1 |  |
| BIRC3 | auto | ok | 12 | 12 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | BIRC2 |  |
| BTK | auto | ok | 12 | 12 | 1 | 0 | 0 |  |  |  | no | 1 | RPL36A-HNRNPH2 |  |
| CBFB | auto | ok | 12 | 12 | 1 | 0 | 0 |  |  |  | no | 1 | SLC6A14 |  |
| EED | auto | ok | 12 | 12 | 2 | 0 | 0 |  |  |  | no | 0.5 | DCLK1 |  |
| IRS1 | auto | ok | 12 | 13 | 2 | 0 | 0 |  |  |  | no | 0.5 | ABCA12 |  |
| RAB35 | auto | ok | 12 | 12 | 3 | 0 | 0 |  |  |  | no | 0.3333 | CDH12 |  |
| SMAD2 | auto | ok | 12 | 13 | 1 | 0 | 0 |  |  |  | no | 1 | RERE |  |
| VEGFA | auto | ok | 12 | 12 | 4 | 0 | 0 |  |  |  | no | 0.5 | ABCC4 |  |
| ALOX12B | auto | ok | 11 | 11 | 2 | 0 | 50 |  |  |  | no | 0.5 | HES7 |  |
| BABAM1 | auto | ok | 11 | 12 | 4 | 0 | 0 |  |  |  | no | 0.5 | ATP6V1E1 |  |
| BCL2 | auto | ok | 11 | 11 | 1 | 0 | 0 |  |  |  | no | 1 | VPS4B |  |
| EPHA7 | auto | ok | 11 | 11 | 2 | 0 | 100 |  |  |  | no | 0.5 | CDH15 |  |
| FANCC | auto | ok | 11 | 11 | 3 | 0 | 0 |  |  |  | no | 0.6667 | CNTLN |  |
| FBXW7 | auto | ok | 11 | 11 | 3 | 0 | 0 |  |  |  | no | 0.3333 | PKD1L1 |  |
| IDH2 | auto | ok | 11 | 11 | 2 | 0 | 0 |  |  |  | no | 0.5 | CHN1 |  |
| SRSF2 | auto | ok | 11 | 11 | 1 | 0 | 100 |  |  |  | no | 1 | NPAS1 |  |
| TNFAIP3 | auto | ok | 11 | 12 | 3 | 0 | 0 |  |  |  | no | 0.3333 | AIG1 |  |
| AURKA | auto | ok | 10 | 10 | 3 | 0 | 66.67 |  |  |  | no | 0.6667 | RAB3GAP2 |  |
| CCND2 | auto | ok | 10 | 10 | 5 | 0 | 60 |  |  |  | no | 0.4 | ETV6 |  |
| EPHB1 | auto | ok | 10 | 10 | 2 | 0 | 50 |  |  |  | no | 0.5 | GMPS |  |
| H3C2 | unresolved | failed | 10 | 10 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| MAPK3 | auto | ok | 10 | 10 | 1 | 0 | 100 |  |  |  | no | 1 | GDPD3 |  |
| PIK3CG | auto | ok | 10 | 10 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| PMS1 | auto | ok | 10 | 11 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| PTPN11 | auto | ok | 10 | 10 | 3 | 0 | 0 |  |  |  | no | 0.3333 | ANKS1B |  |
| SH2B3 | auto | ok | 10 | 10 | 4 | 0 | 25 |  |  |  | no | 0.25 | ATXN2 |  |
| SOX17 | auto | ok | 10 | 10 | 3 | 0 | 0 |  |  |  | no | 0.6667 | XKR4 |  |
| SYK | auto | ok | 10 | 10 | 1 | 0 | 0 |  |  |  | no | 1 | SMC2 |  |
| CHEK1 | auto | ok | 9 | 9 | 2 | 0 | 0 |  |  |  | no | 0.5 | EI24 |  |
| GATA1 | auto | ok | 9 | 9 | 1 | 0 | 0 |  |  |  | no | 1 | HDAC6 |  |
| IFNGR1 | auto | ok | 9 | 9 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | ALDH2 |  |
| MAPKAP1 | auto | ok | 9 | 10 | 1 | 0 | 100 |  |  |  | no | 1 | KIF2A |  |
| PPP6C | auto | ok | 9 | 9 | 1 | 0 | 100 |  |  |  | no | 1 | LCN8 |  |
| PRDM14 | auto | ok | 9 | 9 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| REL | auto | ok | 9 | 9 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| SHQ1 | auto | ok | 9 | 9 | 2 | 0 | 50 |  |  |  | no | 0.5 | FHIT |  |
| AKAP8 | auto | ok | 8 | 8 | 1 | 0 | 100 |  |  |  | no | 1 | NOTCH3 |  |
| AURKB | auto | ok | 8 | 8 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | ALOXE3 |  |
| BARD1 | auto | ok | 8 | 8 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| BCL2L14 | auto | ok | 8 | 8 | 4 | 50 | 0 |  |  |  | no | 0.75 | ETV6 |  |
| CCND3 | auto | ok | 8 | 8 | 2 | 0 | 0 |  |  |  | no | 1 | MED20 |  |
| CD79B | auto | ok | 8 | 8 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | SCN4A |  |
| ERCC3 | auto | ok | 8 | 8 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | CRACDL |  |
| FGF3 | auto | ok | 8 | 8 | 3 | 0 | 0 |  |  |  | no | 0.3333 | PPP6R3 |  |
| GNA11 | auto | ok | 8 | 8 | 1 | 0 | 0 |  |  |  | no | 1 | LMNB2 |  |
| ICOSLG | auto | ok | 8 | 8 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| LYN | auto | ok | 8 | 8 | 4 | 0 | 0 |  |  |  | no | 0.25 | ASIC2 |  |
| NBPF20 | unresolved | failed | 8 | 8 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| PDCD1 | auto | ok | 8 | 8 | 1 | 0 | 0 |  |  |  | no | 1 | LINC01237 |  |
| RRAS | auto | ok | 8 | 8 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | BCL2L12 |  |
| XRCC2 | auto | ok | 8 | 8 | 1 | 0 | 0 |  |  |  | no | 1 | SHH |  |
| ACP3 | auto | ok | 7 | 7 | 7 | 71.43 | 0 |  |  |  | no | 0.4286 | FGFR1 |  |
| AP1B1 | auto | ok | 7 | 7 | 2 | 0 | 0 |  |  |  | no | 1 | EWSR1 |  |
| ASIC2 | auto | ok | 7 | 7 | 6 | 16.67 | 0 |  |  |  | no | 0.3333 | NF1 |  |
| DCBLD1 | auto | ok | 7 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| EPCAM | auto | ok | 7 | 7 | 2 | 0 | 0 |  |  |  | no | 0.5 | GMCL1 |  |
| H3C4 | auto | ok | 7 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| INHA | auto | ok | 7 | 7 | 1 | 0 | 0 |  |  |  | no | 1 | EIF2AK2 |  |
| JUN | auto | ok | 7 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| KLF4 | auto | ok | 7 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| NKX3-1 | auto | ok | 7 | 7 | 1 | 0 | 0 |  |  |  | no | 1 | WRN |  |
| PDPK1 | auto | ok | 7 | 7 | 1 | 0 | 100 |  |  |  | no | 1 | CPPED1 |  |
| RAB11FIP4 | unresolved | failed | 7 | 7 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| RYBP | auto | ok | 7 | 8 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| SDHC | auto | ok | 7 | 7 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | ARHGEF11 |  |
| SHOC2 | auto | ok | 7 | 8 | 1 | 0 | 0 |  |  |  | no | 1 | TSC22D1 |  |
| TENT5C | auto | ok | 7 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| AGAP3 | unresolved | failed | 6 | 6 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| CD274 | auto | ok | 6 | 6 | 1 | 0 | 0 |  |  |  | no | 1 | EYA1 |  |
| CD276 | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| CDH3 | unresolved | failed | 6 | 6 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| FGF4 | auto | ok | 6 | 6 | 2 | 0 | 50 |  |  |  | no | 0.5 | CPT1A |  |
| H1-2 | auto | ok | 6 | 6 | 1 | 0 | 0 |  |  |  | no | 1 | SLC17A3 |  |
| H3-3B | auto | ok | 6 | 6 | 4 | 0 | 0 |  |  |  | no | 0.5 | UNK |  |
| H3C6 | unresolved | failed | 6 | 6 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| HOXB13 | auto | ok | 6 | 6 | 2 | 0 | 0 |  |  |  | no | 0.5 | ABCD3 |  |
| IL7R | auto | ok | 6 | 6 | 2 | 0 | 100 |  |  |  | no | 0.5 | CDH12 |  |
| INHBA | auto | ok | 6 | 6 | 2 | 0 | 50 |  |  |  | no | 0.5 | EPC1 |  |
| ITM2B | auto | ok | 6 | 6 | 3 | 0 | 0 |  |  |  | no | 1 | RB1 |  |
| KANSL1 | auto | ok | 6 | 6 | 5 | 0 | 0 |  |  |  | no | 0.4 | MAP3K14 |  |
| KAZN | auto | ok | 6 | 6 | 4 | 0 | 75 |  |  |  | no | 0.25 | FLT4 |  |
| KLF5 | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| MAX | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| NKX2-1 | auto | ok | 6 | 6 | 2 | 0 | 100 |  |  |  | no | 1 | NKX2-8 |  |
| NOL4 | auto | ok | 6 | 6 | 6 | 66.67 | 0 |  |  |  | no | 0.6667 | FGFR2 |  |
| PDCD1LG2 | auto | ok | 6 | 6 | 2 | 0 | 0 |  |  |  | no | 0.5 | GLIS3 |  |
| PLCD3 | unresolved | failed | 6 | 6 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |
| PMS2 | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| PRDM1 | auto | ok | 6 | 6 | 1 | 0 | 0 |  |  |  | no | 1 | PBK |  |
| PTP4A1 | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| RHOA | auto | ok | 6 | 6 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | RASA2 |  |
| SESN2 | auto | ok | 6 | 7 | 1 | 0 | 0 |  |  |  | no | 1 | TEX46 |  |
| SOCS1 | auto | ok | 6 | 6 | 1 | 0 | 100 |  |  |  | no | 1 | RBM44 |  |
| SPG11 | auto | ok | 6 | 6 | 3 | 0 | 0 |  |  |  | no | 0.6667 | TP53 |  |
| SPRED1 | auto | ok | 6 | 7 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| TSHR | auto | ok | 6 | 6 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| WWTR1 | auto | ok | 6 | 6 | 1 | 0 | 100 |  |  |  | no | 1 | NPM1 |  |
| ZNF598 | auto | ok | 6 | 6 | 5 | 0 | 0 |  |  |  | no | 0.8 | TSC2 |  |
| ZRSR2 | auto | ok | 6 | 7 | 1 | 0 | 0 |  |  |  | no | 1 | SHROOM2 |  |
| AHCYL2 | auto | ok | 5 | 5 | 3 | 0 | 66.67 |  |  |  | no | 1 | SMO |  |
| CASZ1 | auto | ok | 5 | 5 | 5 | 0 | 0 |  |  |  | no | 0.2 | DNMT1 |  |
| EYA2 | auto | ok | 5 | 5 | 4 | 0 | 0 |  |  |  | no | 0.75 | NCOA3 |  |
| H2BC5 | auto | ok | 5 | 5 | 2 | 0 | 0 |  |  |  | no | 0.5 | H1-4 |  |
| H3-5 | auto | ok | 5 | 5 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| KCNU1 | auto | ok | 5 | 5 | 1 | 0 | 0 |  |  |  | no | 1 | RRAS |  |
| KIAA1217 | auto | ok | 5 | 5 | 5 | 80 | 0 |  |  |  | no | 0.8 | FGFR2 |  |
| MED13 | auto | ok | 5 | 5 | 3 | 0 | 66.67 |  |  |  | no | 0.3333 | BRCA1 |  |
| MIPOL1 | auto | ok | 5 | 5 | 3 | 0 | 0 |  |  |  | no | 0.6667 | FOXA1 |  |
| MST1 | auto | ok | 5 | 5 | 1 | 0 | 100 |  |  |  | no | 1 | MON1A |  |
| MUTYH | auto | ok | 5 | 5 | 3 | 0 | 100 |  |  |  | no | 1 | OSBPL9 |  |
| MYOD1 | auto | ok | 5 | 6 | 3 | 0 | 33.33 |  |  |  | no | 0.3333 | DNAJC24 |  |
| NADK | auto | ok | 5 | 5 | 1 | 0 | 0 |  |  |  | no | 1 | CCNL2 |  |
| PPP4R2 | auto | ok | 5 | 5 | 1 | 0 | 0 |  |  |  | no | 1 | RAF1 |  |
| RHEB | auto | ok | 5 | 5 | 2 | 0 | 0 |  |  |  | no | 0.5 | ACTR3B |  |
| SLFN11 | auto | ok | 5 | 5 | 0 | 0 | 0 |  |  |  | no |  |  |  |
| VHL | auto | ok | 5 | 5 | 1 | 0 | 0 |  |  |  | no | 1 | IRAK2 |  |
| WHR1 | auto | ok | 5 | 5 | 2 | 0 | 0 |  |  |  | no | 0.5 | C4A |  |
| XIAP | auto | ok | 5 | 5 | 1 | 0 | 0 |  |  |  | no | 1 | CHRDL1 |  |
| ZFTA | unresolved | failed | 5 | 5 |  |  |  |  |  |  | no |  |  | No canonical transcript/protein could be resolved for this gene. |

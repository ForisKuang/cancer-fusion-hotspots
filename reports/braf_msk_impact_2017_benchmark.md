# BRAF apples-to-apples benchmark: `msk_impact_2017`

Live data were retrieved from the public cBioPortal and Genome Nexus services on
2026-09-03. The [cBioPortal study list](https://www.cbioportal.org/api/studies?projection=SUMMARY&pageSize=10000000&pageNumber=0)
identifies `msk_impact_2017` as **MSK-IMPACT
Clinical Sequencing Cohort (MSK, Nat Med 2017)** with 10,945 samples. Its molecular
profile list for [`msk_impact_2017`](https://www.cbioportal.org/api/studies/msk_impact_2017/molecular-profiles)
includes `msk_impact_2017_structural_variants` (alteration type
`STRUCTURAL_VARIANT`, datatype `SV`), confirming that structural-variant data are
present.

## Results

- BRAF structural variants returned: **48**
- Resolved protein fusions: **35/35** (no protein-fusion record was skipped)
- In-frame protein fusions: **33/35 (94.3%)**
- PF07714 kinase-domain retained overall: **31/35 (88.6%)**
- PF07714 kinase-domain retained among in-frame fusions: **31/33 (93.9%)**
- In-frame and PF07714-retained jointly: **31/35 (88.6%)**
- Unique partner genes: **19**
- One-sided Fisher exact test: **odds ratio = infinity, p = 0.0100840**
- Frame/domain table `[[retained/in-frame, retained/other], [not-retained/in-frame,
  not-retained/other]]`: **`[[31, 0], [2, 2]]`**
- Breakpoint-permutation empirical p-value (10,000 seeded permutations): **0.328767**

## Side-by-side comparison

| Metric | Zehir et al. (PMC5461196) finding | `msk_impact_2017` live reanalysis | `msk_impact_50k_2026` committed benchmark |
|---|---:|---:|---:|
| BRAF SVs returned | Not reported by the cited 33-fusion finding | 48 | 251 |
| Resolved BRAF protein fusions | 33 | 35 | 174 |
| In-frame | 33/33 (100%) | 33/35 (94.3%) | 151/174 (86.8%) |
| PF07714 retained overall | 33/33 (100%) | 31/35 (88.6%) | 161/174 (92.5%) |
| PF07714 retained among in-frame | 33/33 (100%) | 31/33 (93.9%) | 142/151 (94.0%) |
| In-frame and PF07714-retained jointly | 33/33 (100%) | 31/35 (88.6%) | 142/174 (81.6%) |
| Unique partner genes | Not reported by the cited finding | 19 | 63 |
| Fisher odds ratio (one-sided) | Not applicable to the reported all-positive subset | infinity | 3.32164 |
| Fisher p-value (one-sided) | Not applicable to the reported all-positive subset | 0.0100840 | 0.0737799 |

## Interpretation

The older cohort is closer to the paper on the in-frame rate and on the joint
in-frame-plus-retained endpoint than the 50k successor cohort. It nevertheless does
**not** reproduce the paper's clean 33/33 result: the live portal currently returns
35 BRAF protein-fusion records, including two with unknown frame, and four are
classified as PF07714-lost (two of those four are in-frame). Overall domain retention
is actually lower than in the 50k benchmark (88.6% versus 92.5%), while domain
retention among in-frame fusions is essentially the same (93.9% versus 94.0%). Thus
the 2017 cohort is a closer reproduction only for frame/joint purity, not uniformly
across every retention metric, and neither live analysis yields 100%/100%.

The difference between the paper's 33 reported fusions and the 35 live cBioPortal
protein-fusion records is reported as observed; no records were removed to force the
published denominator.

The reported infinite odds ratio is a small-sample zero-cell artifact, not a modeling error: the frame/domain contingency table (`[[31, 0], [2, 2]]`) has no observed out-of-frame-or-unknown-frame event that also retains the kinase domain, so the odds ratio's denominator is zero. Fisher's exact p-value remains well-defined and is reported directly from the hypergeometric distribution regardless of the odds-ratio degeneracy.

## Method and provenance

The existing shared CLI was run without a study-specific code fork:

```bash
cfh real-benchmark BRAF msk_impact_2017 \
  --output-dir /tmp/cfh-msk2017-live \
  --n-permutations 10000
```

The pipeline queried cBioPortal's
`msk_impact_2017_structural_variants` profile using the configured BRAF Entrez gene
ID (673), selected records explicitly annotated as `Protein Fusion`, adapted them to
the production SV schema, and normalized frame and 5'/3' orientation from the source
annotations. BRAF genomic breakpoints were then mapped to the Genome Nexus canonical
transcript and evaluated with the existing domain-retention algorithm against
PF07714 (amino acids 458-712). Counts are event-level with no patient deduplication;
the Fisher `other` column combines out-of-frame and unknown-frame events as defined
by the existing algorithm.

The study and profile identities were independently verified through the cBioPortal
REST study-list and study molecular-profile endpoints before running the benchmark.

## Partners

AGAP3 (2), AGK (3), CCDC6 (1), CDK5RAP2 (2), CUL1 (1), FAM131B (1), GIPC2 (1),
KIAA1549 (4), MKRN1 (3), OSBPL9 (1), PARP12 (1), PHTF2 (1), PJA2 (1), PRKAR1B (1),
PRKAR2B (1), RBM33 (1), SCRIB (1), SND1 (8), ZNF207 (1)

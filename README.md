# Analysis code for: A peroxisome-related prognostic signature is associated with reduced lymphocytic infiltration in hepatocellular carcinoma

Reproducible analysis code for the manuscript submitted to *PLOS ONE* (2026).

**v1.1** adds the analyses performed in response to peer review; see

[CHANGELOG.md](CHANGELOG.md) for the exact list.

## Overview

- Single-cell RNA-seq processing of five HCC datasets (**76,762 cells** after QC)
- A **131-gene peroxisome-related gene (PRG)** set used to score peroxisome

  activity (PXS), and a **ten-gene peroxisome-related prognostic signature (PRS)**

  derived within an eight-algorithm framework
- Apparent TCGA-LIHC C-index **0.761**, which falls to **0.556-0.564** under

  leakage-free nested cross-validation
- External evaluation in four GEO cohorts and ICGC-LIRI-JP; random-effects

  meta-analysis pooled HR **1.57** (95% CI 1.18-2.09), I-squared **0%**
- Multivariable Cox regression adjusted for age, sex and AJCC stage: PRS

  HR **2.66 per SD** (95% CI 2.06-3.44), C-index **0.61 -> 0.78**
- Immune, mutational and drug-sensitivity associations, with

  Benjamini-Hochberg correction and a tumour-purity sensitivity analysis
- Figure generation (main Figs 1-8) and supporting figures/tables

We report the negative and the positive results with equal weight: the

single-cell difference between malignant and non-malignant hepatocytes does not

survive patient-level analysis (58 paired patients, p = 0.54), and the

ten-gene score performs at the level of the null under nested cross-validation.

All data used are publicly available; see **Data availability**.

## Repository structure

```
.
|-- README.md
|-- CHANGELOG.md
|-- LICENSE                  # MIT
|-- requirements.txt         # pinned dependency versions
|-- .gitignore
|-- config.py                # paths, random seed, cohort accessions
|-- scripts/                 # flat, so that intra-package imports resolve
|   |-- figstyle.py          # shared plotting style
|   |-- pubs.py              # citation helper
|   |-- core.py              # shared utilities + ML framework
|   |-- fig_bulk.py          # Figs 2-5, 8
|   |-- fig_sc.py            # Figs 1, 6, 7
|   |-- analysis_tcga.py     # PRS construction, TCGA associations, CV table
|   |-- analysis_bulk2.py    # bulk association and external validation
|   |-- build_prg.py         # 131-gene PRG set and PXS scoring
|   |-- sc_preprocess.py / sc_analysis.py / sc_annotate.py
|   |-- build_supp_figures.py / build_supp_tables.py
|   |-- revision_*.py        # analyses added in the first revision
|   `-- r1_*.py              # analyses added in the second revision
`-- data/
    `-- processed/           # lightweight outputs only (see data/processed/README.md)
```

Scripts are kept flat so that `import figstyle`, `import core` and

`from config import PROJECT_ROOT` resolve without modification. Every hardcoded

development path has been redirected through `config.PROJECT_ROOT`, so the

repository is relocatable.

## Environment

Python **3.13.12**:

```bash
pip install -r requirements.txt
```

## Data availability

No new data were generated. Raw data are public; the ~18 GB of raw files and the

bulk-expression matrices are **not** bundled. Download them from the accessions

below into `data/raw/`, or re-run the acquisition scripts (`build_geo.py`,

`build_tcga.py`, `build_icgc.py`, `fetch_*.py`). MSigDB gene sets require a free

license and go in `data/genesets/`.

| Dataset      | Repository and direct link                                     | Use                                                                   |
| ------------ | -------------------------------------------------------------- | --------------------------------------------------------------------- |
| GSE149614    | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614> | scRNA-seq atlas                                                       |
| GSE151530    | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151530> | scRNA-seq atlas                                                       |
| GSE125449    | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449> | scRNA-seq atlas                                                       |
| GSE146115    | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146115> | scRNA-seq atlas                                                       |
| GSE98638     | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638>  | scRNA-seq atlas                                                       |
| GSE14520     | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520>  | validation with survival                                              |
| GSE76427     | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76427>  | validation with survival                                              |
| GSE116174    | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116174> | validation with survival                                              |
| GSE10186     | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10186>  | validation with survival                                              |
| TCGA-LIHC    | <https://portal.gdc.cancer.gov/projects/TCGA-LIHC>             | training, KM (Firehose 2016_01_28, <https://gdac.broadinstitute.org>) |
| ICGC-LIRI-JP | <https://dcc.icgc.org/projects/LIRI-JP>                        | validation (release 28)                                               |
| GDSC1        | <https://www.cancerrxgene.org/downloads/bulk_download>         | drug sensitivity                                                      |

A lightweight subset of processed outputs is bundled under `data/processed/`

(clinical and risk tables, PRS coefficients, the PRG set, single-cell

annotations, and the revision-round result tables) so that the figure and table

scripts can be re-run without the full matrices.

## How to reproduce

```bash
# 0. fetch raw data (or place it manually under data/raw/)
python scripts/build_geo.py
python scripts/build_tcga.py
python scripts/build_icgc.py

# 1. single-cell pipeline + Figs 1, 6, 7
python scripts/sc_preprocess.py
python scripts/sc_analysis.py
python scripts/sc_annotate.py
python scripts/fig_sc.py

# 2. PRS construction, TCGA associations, ten-fold CV table
python scripts/analysis_tcga.py        # writes data/processed/TCGA_PRS_coefficients.csv

# 3. bulk association, external validation, Figs 2-5, 8
python scripts/analysis_bulk2.py
python scripts/fig_bulk.py

# 4. first-revision re-analyses
python scripts/revision_reanalyses.py      # BH correction, nested CV, age/sex model
python scripts/revision_sc_reanalyses.py   # variance decomposition, patient-level tests
python scripts/revision_bulk_reanalyses2.py # external rigour, meta-analysis, PH test

# 5. second-revision re-analyses
python scripts/r1_multivariable_stage.py   # multivariable Cox with AJCC stage
python scripts/r1_superpc_cv.py            # SuperPC entry for Fig 2B
python scripts/r1_sensitivity_hcc_only.py  # HCC-only GSE151530 (needs the annotated h5ad)

# 6. supporting figures and tables
python scripts/build_supp_figures.py
python scripts/build_supp_tables.py
```

Steps 4-5 need the intermediate matrices produced in steps 0-3; they are the

only steps that require the full download. Runtime is dominated by

single-cell and bulk preprocessing (hours on a multicore node); the

PRS/figure/table steps take minutes.

## What each revision script answers

| Script                         | Output                                                                                                                                               | Purpose                                                                                                                                          |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `revision_reanalyses.py`       | `revision_fdr.csv`, `revision_nestedCV_cindex.csv`, `revision_multivariable_cox.csv`, `revision_multivariable_summary.json`                          | Benjamini-Hochberg correction; leakage-free nested cross-validation; age/sex-adjusted Cox                                                        |
| `revision_sc_reanalyses.py`    | `revision_sc_variance.json`, `revision_sc_batch.csv`, `revision_sc_patient_level.csv`, `revision_sc_patient_tests.csv`                               | variance decomposition of the peroxisome score; patient-level pseudobulk tests                                                                   |
| `revision_bulk_reanalyses2.py` | `revision_external_rigor.csv`, `revision_meta_heterogeneity.json`, `revision_ph_test.csv`, `revision_purity_confounding.csv`, `revision_n_table.csv` | per-cohort C-index with bootstrap CI, random-effects meta-analysis, proportional-hazards test, tumour-purity sensitivity, sample-size accounting |
| `r1_multivariable_stage.py`    | `revision_multivariable_stage.csv/.json`                                                                                                             | nested models M1 (age, sex) -> M2 (+AJCC stage) -> M3 (+PRS), likelihood-ratio test, Child-Pugh sensitivity                                      |
| `r1_superpc_cv.py`             | `TCGA_ML_CV_cindex_v2.csv`, `superpc_cv.json`                                                                                                        | recovers the SuperPC fold with a convergent implementation; original table left intact                                                           |
| `r1_sensitivity_hcc_only.py`   | `revision_sc_hcc_only*.csv/.json`                                                                                                                    | restricts GSE151530 to its 32 HCC samples (14 are iCCA)                                                                                          |

Supporting tables S13 and S14 are assembled directly from the `revision_*` CSV


files listed above; no additional script is required.

## The PRS (ten genes, TCGA-derived coefficients)

| Gene     | Coefficient | Direction  |
| -------- | ----------- | ---------- |
| IDH1     | +0.467      | risk       |
| UBE2D2   | +0.383      | risk       |
| AOC1     | +0.342      | risk       |
| SLC25A17 | +0.204      | risk       |
| UBB      | -0.174      | protective |
| CRAT     | -0.189      | protective |
| MTARC2   | -0.247      | protective |
| ABCD2    | -0.282      | protective |
| ACOXL    | -0.333      | protective |
| PECR     | -0.433      | protective |

Full coefficients: `data/processed/TCGA_PRS_coefficients.csv`; gene set:

`data/processed/peroxisome_gene_set.json`.

## License

MIT (see `LICENSE`). PLOS ONE requires author-generated code to be available

without restriction; MIT satisfies this.

## Citation

> Jin G, Zhao G, You J, Li M. A peroxisome-related prognostic signature is
>
> associated with reduced lymphocytic infiltration in hepatocellular carcinoma.
>
> *PLOS ONE* (2026). DOI: 10.1371/journal.pone.XXXXXXX

Code archive (all versions): <https://doi.org/10.5281/zenodo.22640219>

Version v1.1: <https://doi.org/10.5281/zenodo.NNNNNNNN>

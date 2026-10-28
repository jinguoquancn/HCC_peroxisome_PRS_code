# data/processed - lightweight subset

Only the lightweight outputs needed to inspect the results and to run the
figure/table scripts are included. Raw data and bulk-expression matrices are
**not** bundled: regenerate them with the pipeline or download them from the
public repositories listed in the top-level README.

## Bundled here
- `peroxisome_gene_set.json` - 131-gene peroxisome-related gene (PRG) set
- `TCGA_PRS_coefficients.csv` - ten-gene PRS coefficients
- `*_clinical.csv`, `*_validation_surv.csv`, `*_risk.csv` - cohort tables
- `TCGA_*` - TCGA-LIHC derived tables (mutations, TMB, DEG, PXS, cross-validation
  C-index, clinical with PRS/PXS)
- `sc/` - single-cell annotations, CNV, cell-cell communication, pseudotime
- association outputs: `immune_infiltration_PRS.csv`, `mutation_PRS.csv`,
  `drug_sensitivity_PRS.csv`, `signature_comparison.csv`, `GSEA_PRS_high_vs_low.csv`
- `revision_*` - result tables added in revision rounds (see CHANGELOG.md)


# TAD–PDAC analysis scripts

R and Python scripts accompanying the manuscript on new-onset or worsening diabetes, survival, and intratumoral T-cell densities in resected pancreatic ductal adenocarcinoma (PDAC).

The analyses include 162 patients in the surgical cohort (64 TAD, 98 non-TAD), 36 patient-derived organoids (6 TAD, 30 non-TAD), and 163 TCGA-PAAD primary tumors. TAD denotes new-onset or worsening diabetes within three years before PDAC diagnosis.

## Contents

| Script | Analysis | Manuscript outputs |
| --- | --- | --- |
| `01_survival_analysis.R` | Cox models, survival comparisons, sensitivity and interaction analyses | Table 2; Fig. 1b–c; Supplementary Figs. S2–S3; Tables S4–S6A and survival rows of S6B |
| `02_deseq2_differential_expression.R` | Organoid differential expression, PCA, platform sensitivity | Fig. 3a, 3c–d; Supplementary Fig. S5; Table S9A |
| `03_preranked_gsea.R` | Hallmark and Reactome enrichment | Fig. 3b; Supplementary Tables S8A–S8B |
| `04_secreted_factor_tpm.py` | TPM-based sensitivity analyses | Supplementary Table S9B |
| `05_tcell_adjusted_models.R` | Adjusted T-cell density models and sensitivity analyses | Fig. 2c; T-cell rows of Supplementary Table S6B; Tables S6C, S7B–S7D |
| `06_tcga_paad_analysis.R` | GDF15-related expression, tumor-purity adjustment, survival | Fig. 4c–e; Supplementary Table S11; Supplementary Fig. S6 |
| `07_driver_gene_wgs_analysis.R` | Organoid WGS-derived driver alteration frequencies | Supplementary Fig. S4 |

Scripts generate numerical results and figure source data rather than final figure graphics.

## Requirements

- R with `readxl`, `dplyr`, `tibble`, `openxlsx`, `survival`, `broom`, `sandwich`, `msigdbr`, `DESeq2`, `apeglm`, and `clusterProfiler` (including `DOSE` and `fgsea`).
- Python with `pandas`, `numpy`, `scipy`, and `openpyxl`.
- MSigDB release **2025.1.Hs** for GSEA. Exact output may depend on package versions.

## Input data

Run scripts from the repository root. File structures, column coding, and sample requirements are described in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md).

| Local path | Scripts |
| --- | --- |
| `data/clinical_data.xlsx` | 01, 05 |
| `data/count_matrix.xlsx` | 02 |
| `data/organoid_expression_TPM.xlsx` | 04 |
| `data/organoid_driver_alterations.xlsx` | 07 |
| `data/TCGA/TCGA-PAAD.star_fpkm-uq.tsv.gz` | 06 |
| `data/TCGA/TCGA-PAAD.clinical.tsv.gz` | 06 |
| `data/TCGA/TCGA-PAAD.survival.tsv.gz` | 06 |
| `data/TCGA/gencode.v36.annotation.gtf.gene.probemap` | 06 |
| `data/TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | 06 |

TCGA-PAAD data were obtained from the UCSC Xena GDC Hub and PanCanAtlas/Genomic Data Commons resources.

## Execution

```bash
Rscript 01_survival_analysis.R
Rscript 02_deseq2_differential_expression.R
Rscript 03_preranked_gsea.R
python 04_secreted_factor_tpm.py
Rscript 05_tcell_adjusted_models.R
Rscript 06_tcga_paad_analysis.R
Rscript 07_driver_gene_wgs_analysis.R
```

Run `02` before `03` and `04`: both downstream scripts use DESeq2 output. Script `07` optionally checks organoid sample IDs against metadata produced by `02`. Other scripts can be run independently with the corresponding data in place.

Results are written under `results/`; software versions are recorded in the generated logs or manifests. The 13 Table 2 univariable analyses require the validated input fields listed in the data dictionary.

## Analysis conventions

- **Survival:** the primary Cox model includes TAD, pathologic stage, resection margin, categorized NLR, and categorized CA19-9. The OS sensitivity analyses also include log2-transformed intratumoral CD8 density.
- **T-cell densities:** log-linear regression with HC3 robust inference and BH adjustment across the three primary markers (CD4, CD8, FOXP3).
- **Organoid RNA-seq:** DESeq2 primary design `~ DM_group`; apeglm-shrunken log2 fold changes with unshrunken Wald-test P values.
- **GSEA:** `sign(shrunken log2FC) × −log10(nominal P)`; MSigDB 2025.1.Hs; BH-adjusted P <0.25 for Hallmark (S8A) and <0.05 for Reactome (S8B).
- **TPM:** two-sided Mann–Whitney U tests without multiplicity adjustment; supplementary sensitivity results for Fig. 3d genes.
- **Organoid driver genes:** two-sided Fisher exact tests of pre-derived alteration calls; no raw WGS processing.
- **TCGA-PAAD:** GDF15-related transcriptomic analyses and ABSOLUTE purity sensitivity. The supplied expression matrix is already log2-transformed and must not be transformed again.

## Data availability

Clinical patient-level data are not publicly distributed for privacy and ethics reasons. Organoid RNA-seq and whole-genome sequencing data are available under controlled access at the European Genome-phenome Archive (EGA): `EGAS00001007212` and `EGAS00001007211`, respectively.

Restricted-data analyses require authorized access to the corresponding inputs. Stage IV-adjusted DESeq2/GSEA analyses (Supplementary Table S8C) are not included because donor-level stage annotations are unavailable for public release.

# TAD-PDAC analysis scripts

Analysis scripts accompanying the manuscript:

**New-onset or worsening diabetes, survival, and intratumoral T-cell densities in resected pancreatic ductal adenocarcinoma**

## Scope

This repository contains the statistical and bioinformatic analysis scripts used for the manuscript and its Supplementary Information.

Patient-level clinical data are not distributed because of participant privacy and ethical restrictions. Patient-derived organoid sequencing data are available through controlled access at the European Genome-phenome Archive:

- RNA sequencing: **EGAS00001007212**
- Whole-genome sequencing: **EGAS00001007211**

Public TCGA-PAAD transcriptomic and clinical data were obtained from the UCSC Xena GDC Hub. PanCanAtlas ABSOLUTE tumour-purity estimates were obtained from the Genomic Data Commons.

Clinicopathological comparisons and selected descriptive/immunohistochemical group comparisons were performed in **IBM SPSS Statistics version 29.0.1**. SPSS syntax is not included in this repository. The adjusted T-cell models and the other statistical/bioinformatic analyses described below are represented by executable R or Python scripts.

## Repository structure

```text
TAD-PDAC-analysis-scripts/
├── README.md
├── DATA_DICTIONARY.md
├── .gitignore
├── 01_survival_analysis.R
├── 02_deseq2_differential_expression.R
├── 03_preranked_gsea.R
├── 04_secreted_factor_tpm.py
├── 05_tcell_adjusted_models.R
├── 06_tcga_paad_analysis.R
└── 07_driver_gene_wgs_analysis.R
```

The patient-level clinical input file (`data/clinical_data.xlsx`) is intentionally not distributed.

## Active scripts

| Script | Analysis | Manuscript mapping |
|---|---|---|
| `01_survival_analysis.R` | OS/CSS, primary Cox models, sensitivity models, four diabetes phenotypes, diabetes-only analyses, stage-stratified models, and TAD-by-stage interaction | Fig. 1; Table 2; Supplementary Figs. S1–S2; Supplementary Tables S3–S5 |
| `02_deseq2_differential_expression.R` | Organoid DESeq2 differential expression, PCA, and sequencing-platform sensitivity | Fig. 3a, c, d; Supplementary Fig. S4; Supplementary Table S8A |
| `03_preranked_gsea.R` | Hallmark and Reactome pre-ranked GSEA | Fig. 3b; Supplementary Table S7 |
| `04_secreted_factor_tpm.py` | TPM-based sensitivity analysis of selected organoid genes | Supplementary Table S8B |
| `05_tcell_adjusted_models.R` | Adjusted CD4/CD8/FOXP3 models, HC3 inference, bootstrap CIs, HbA1c and body-status sensitivities, four-phenotype analyses, and stage-interaction analyses | Fig. 2c; Supplementary Tables S5B–C and S6B–D |
| `06_tcga_paad_analysis.R` | TCGA-PAAD GDF15/immune/CYT/DDIT3 analyses, purity adjustment, histology-restricted sensitivity, and OS | Fig. 4c–e; Supplementary Fig. S5; Supplementary Table S10 |
| `07_driver_gene_wgs_analysis.R` | WGS-derived KRAS/TP53/CDKN2A/SMAD4 comparisons in organoids | Supplementary Fig. S3 |

## Key reproducibility decisions

### Survival analyses

The primary multivariable survival model includes:

`TAD + pathologic stage + resection margin status + dichotomised NLR + dichotomised CA19-9`

with NLR dichotomised at the cohort median (2.1) and CA19-9 dichotomised at 37 U/mL.

Sensitivity analyses include:

- additional adjustment for continuous preoperative HbA1c;
- replacement of TAD with continuous HbA1c;
- additional adjustment for adjuvant chemotherapy;
- additional adjustment for continuous age;
- additional adjustment for continuous surgery year;
- replacement of dichotomised CA19-9 with `log2(CA19-9)`;
- `log2(CA19-9)` plus adjuvant chemotherapy;
- exclusion of patients who received neoadjuvant chemotherapy;
- classification-restricted analyses excluding treatment-only worsening diabetes;
- diabetes-only classification sensitivities;
- additional adjustment for BMI, SMI, serum albumin and ECOG-PS.

A further **OS-only** sensitivity model treats both NLR and CA19-9 as continuous variables:

`TAD + pathologic stage + resection margin status + continuous NLR + log2(CA19-9)`

The TAD estimate from this model is:

**HR 1.93 (95% CI 1.31–2.85), P = 0.001.**

This model corresponds to **Supplementary Table S4 (overall survival)** and is **not included in Supplementary Table S3 (cancer-specific survival)**.

To reduce the risk of mixing endpoints, `01_survival_analysis.R` writes separate output files:

- `supplementary_table_S3_css_sensitivity.csv`
- `supplementary_table_S4_os_sensitivity.csv`

The TAD-by-stage interaction P values are taken from the **Wald test for the interaction coefficient**. The script also checks the stage-specific sample sizes and event counts reported in Supplementary Fig. S2.

### T-cell models

The main adjusted T-cell models use log-transformed cell density and adjust for:

`TAD + age + pathologic stage + log2(CA19-9) + continuous NLR + neoadjuvant chemotherapy`

HC3 heteroscedasticity-robust standard errors are used. Adjusted geometric-mean ratios with 95% confidence intervals and nominal P values are reported for CD4+, CD8+ and FOXP3+ densities; Benjamini–Hochberg-adjusted q values across the three markers are additionally reported as a multiplicity assessment.

Percentile bootstrap confidence intervals use 2,000 patient-level resamples. Supporting four-category T-cell analyses use patients without diabetes as the reference group.

### Organoid RNA-seq

The primary DESeq2 design is `~ DM_group`, with non-TAD as the reference group. Genes with row-summed raw counts `<=10` are excluded. Log2 fold changes are shrunken using `apeglm`. Differential expression is defined by BH-adjusted `P < 0.05`, with no fold-change threshold.

The sequencing-platform sensitivity model uses:

`~ platform + DM_group`

PCA uses variance-stabilised counts.

### GSEA

Pre-ranked GSEA uses:

`sign(apeglm-shrunken log2FC) * -log10(nominal DESeq2 P value)`

Ties are ordered deterministically. Hallmark and Reactome collections are evaluated with gene-set sizes of 15–500 and random seed 42. The prespecified enrichment threshold is `q < 0.25`.

### TPM sensitivity analyses

TPM-based comparisons use `log2(TPM + 1)` and two-sided Mann–Whitney U tests. Because `log2(TPM + 1)` is monotonic, the transformation does not change the rank ordering.

### TCGA-PAAD

The TCGA workflow:

- restricts samples to primary-tumour `-01A` aliquots;
- excludes `_PAR_Y` records before removing Ensembl version suffixes;
- maps genes using GENCODE v36;
- uses the Xena FPKM-UQ expression scale as supplied by the selected input matrix, without an additional transformation;
- evaluates GDF15 as median-dichotomised and continuous;
- evaluates CD8A, GZMA, GZMB, PRF1, CYT and DDIT3;
- applies BH correction to the prespecified outcome families;
- adjusts partial Spearman correlations for ABSOLUTE tumour purity;
- uses Fisher-z confidence intervals with standard error `1 / sqrt(n - k - 3)`;
- includes a histology-restricted sensitivity analysis.

TCGA-PAAD lacks the longitudinal glycaemic and treatment information required to define TAD and is therefore used as an independent transcriptomic reference rather than as external validation of the TAD classification.

## Software

The manuscript reports:

- R 4.3.3
- Python 3.10.12
- IBM SPSS Statistics 29.0.1

Individual scripts record package/session information where applicable.

## Code availability

Analysis scripts used in this study are publicly available in this repository. Patient-level clinical data are not included because of participant privacy and ethical restrictions.

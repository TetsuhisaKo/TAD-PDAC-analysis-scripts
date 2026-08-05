# TAD-PDAC analysis scripts

Analysis scripts accompanying the manuscript:

**Clinically defined tumor-associated diabetes, intratumoral T-cell densities, and survival in resected pancreatic ductal adenocarcinoma**

## Scope

This repository reproduces the **statistical and bioinformatic results** reported
in the current manuscript. Figures were generated separately and **plotting code
is not included**.

Patient-level clinical data are not distributed because of participant privacy
and ethical restrictions. RNA-sequencing data generated in this study are
available through the European Genome-phenome Archive under accession
**EGAS00001007212**, and whole-genome sequencing data under accession
**EGAS00001007211**. Public TCGA-PAAD data were obtained from the UCSC Xena GDC
Hub.

Clinicopathological comparisons and selected descriptive/IHC group comparisons
were performed in **IBM SPSS Statistics version 29.0.1**. SPSS syntax is not
included in this repository. In particular, the unadjusted Fig. 2a comparisons
and the Fig. 4b GDF15 IHC group comparison are not represented by executable
SPSS code here. The adjusted Fig. 2c T-cell models are R analyses and are
included in `05_tcell_adjusted_models.R`.

## Active scripts

| Script | Analysis |
|---|---|
| `01_survival_analysis.R` | OS/CSS; sensitivity models; four diabetes phenotypes; stage-stratified HRs; Wald TAD-by-stage interaction |
| `02_deseq2_differential_expression.R` | Organoid DESeq2; explicit 500-gene PCA; platform-adjusted sensitivity |
| `03_preranked_gsea.R` | Hallmark/Reactome pre-ranked GSEA; deterministic tie ordering; MSigDB/version logging |
| `04_secreted_factor_tpm.py` | TPM-based sensitivity analysis of selected Fig. 3e genes |
| `05_tcell_adjusted_models.R` | Fig. 2c HC3 models; 2,000-resample bootstrap; HbA1c sensitivity; CD45RO; phenotype and stage-interaction analyses |
| `06_tcga_paad_analysis.R` | TCGA GDF15/immune/CYT/DDIT3; purity adjustment; duct-only sensitivity; OS |
| `07_driver_gene_wgs_analysis.R` | Organoid KRAS/TP53/CDKN2A/SMAD4 Fisher tests |

The files are numbered consecutively so that there are no unexplained gaps.
Analyses removed from the revised manuscript, including the former exploratory
GDF15 ELISA/protein-mRNA and ER-stress leave-one-out workflows, are not part of
this active set.

## Key reproducibility decisions

### Survival interaction

The TAD-by-stage interaction P values reported for OS and CSS are taken from the
**Wald test for the interaction coefficient**. The clinical script also asserts
the stage-specific sample sizes and event counts used in Supplementary Fig. S2.

### T-cell models

The primary immune models use log-transformed cell density and adjust for age,
pathologic stage, log2-transformed CA19-9, continuous NLR, and neoadjuvant
chemotherapy. HC3 heteroscedasticity-robust covariance is used, and
Benjamini-Hochberg correction is applied across the three primary outcomes
(CD4, CD8, FOXP3). Percentile bootstrap confidence intervals use 2,000
patient-level resamples with a fixed random seed.

### PCA

Supplementary Fig. S4 uses variance-stabilized RNA-seq counts and the **500 most
variable genes**. The `ntop=500` setting is explicit in the script.

### GSEA

The GSEA ranking metric is
`sign(shrunken log2FC) * -log10(nominal P)`.

Ties are ordered deterministically. The script records the installed
`clusterProfiler` and `msigdbr` versions and, when exposed by `msigdbr`, the
MSigDB database release.

### TPM Mann-Whitney analyses

The TPM sensitivity analysis uses `log2(TPM+1)`. This is a strictly monotonic
transformation and therefore does not change the Mann-Whitney ranks or P value.
The transformation is retained only for consistency with the expression-scale
description.

### TCGA

The TCGA script:

- explicitly restricts sample barcodes to `-01A`;
- removes `_PAR_Y` records before stripping Ensembl version suffixes;
- records input-file MD5 checksums;
- records the Ensembl IDs actually used;
- uses the already transformed Xena FPKM-UQ matrix directly, with no second log transform;
- defines GDF15-high as expression **at or above the median**;
- applies BH correction across the six prespecified median-split outcomes;
- uses the corrected Fisher-z standard error for partial correlation:
  `1 / sqrt(n - k - 3)`.

## Software

The manuscript reports:

- R 4.3.3
- Python 3.10.12
- IBM SPSS Statistics 29.0.1

Individual scripts write session/package-version information where applicable.

## Code availability

Suggested manuscript wording:

> Analysis scripts used in this study are publicly available at https://github.com/TetsuhisaKo/TAD-PDAC-analysis-scripts. Clinicopathological comparisons and selected unadjusted immunohistochemical group comparisons, including those underlying Fig. 2a and Fig. 4b, were performed using IBM SPSS Statistics version 29.0.1 and are not represented by executable code in the repository. Patient-level clinical data are not included because of participant privacy and ethical restrictions.

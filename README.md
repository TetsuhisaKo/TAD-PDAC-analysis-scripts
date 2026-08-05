Manuscript title:Clinically defined tumor-associated diabetes, intratumoral T-cell densities, and survival in resected pancreatic ductal adenocarcinoma

This bundle revises the scripts supplied for the manuscript so that file descriptions,analysis families, cohort definitions, multiple-testing procedures, and figure/tablemapping match the current manuscript and supplementary material.

Active scripts

01_survival_analysis.R

Overall survival (OS) and cancer-specific survival (CSS)

Primary multivariable model:TAD + pathologic stage + resection margin + NLR + CA19-9

Sensitivity analyses:continuous HbA1c; adjuvant chemotherapy; log2-continuous CA19-9;log2-continuous CA19-9 + adjuvant chemotherapy; exclusion of neoadjuvant cases

Four-category diabetes analysis

Diabetes-only TAD vs long-standing stable diabetes, with and without HbA1c

TAD-by-stage interactions

Stage-stratified adjusted HRs for Supplementary Fig. S2

02_deseq2_differential_expression.R

36 patient-derived PDAC organoids: TAD n=6, non-TAD n=30

Primary DESeq2 design: ~ DM_group

Low-count filter: row-summed raw count >10

apeglm-shrunken log2 fold changes

BH-adjusted P<0.05, no fold-change threshold

Platform-adjusted sensitivity analysis: ~ platform + DM_group

Variance-stabilized PCA coordinates

Selected-gene output for Supplementary Table S7A

03_preranked_gsea.R

Pre-ranked GSEA using sign(log2FC) * -log10(P)

Hallmark and Reactome collections

Gene-set size 15–500; random seed 42

Prespecified GSEA threshold q<0.25

Supplementary Table S6:

Hallmark panel: all nominal P<0.05

Reactome panel: all q<0.25

04_secreted_factor_tpm.py

TPM-based sensitivity analysis for the 12 genes shown in Fig. 3e

Test is performed on log2(TPM+1) using a two-sided asymptoticMann–Whitney U test

Descriptive medians and IQRs are reported on the original TPM scale

Produces the current Supplementary Table S7B layout

GDF15 expected P≈0.040; IL33 expected P≈0.155

05_tcga_gdf15_analysis.R

Histologically restricted TCGA-PAAD expression cohort: n=163

infiltrating duct carcinoma, NOS: n=143

adenocarcinoma, NOS: n=20

OS-evaluable cohort: n=162

GDF15 median split is defined in n=163 before OS restriction

Six median-split comparisons:CD8A, GZMA, GZMB, PRF1, CYT, DDIT3

BH correction across those six comparisons

Continuous Spearman correlations for the same six outcomes

ABSOLUTE tumor-purity-adjusted partial Spearman correlations:n=143, with BH correction across six

Histologically stringent duct-carcinoma sensitivity analysis:n=143; purity-evaluable n=123

TCGA OS: GDF15-low n=81, high n=81; log-rank P≈0.75

Maps to Fig. 4c–e, Supplementary Table S8, Supplementary Fig. S5

08_driver_gene_wgs_analysis.R

Fisher exact tests for WGS-derived organoid alteration calls

KRAS, TP53, CDKN2A, SMAD4

TAD n=6, non-TAD n=30

Maps to Supplementary Fig. S3

Removed from the active final analysis

09_er_upr_leave_one_out_analysis.R is not part of the current manuscript.The current Supplementary Table S7 contains:

panel A: primary vs sequencing-platform-adjusted DESeq2 results

panel B: TPM-based sensitivity analysis

The previous leave-one-out ER/UPR overlap analysis should therefore not be describedas Supplementary Table S7 in the final GitHub repository. A note is retained inarchive/09_er_upr_leave_one_out_analysis_ARCHIVED.txt.

Important TCGA expression-scale note

The UCSC Xena TCGA-PAAD.star_fpkm-uq.tsv.gz matrix used for the final analysescontains the log2-transformed FPKM-UQ values used in the manuscript. The analysisscript therefore uses those values directly and does not apply a secondlog2(x+1) transformation.

Controlled-access data

Patient-level clinical and organoid data are not included in this bundle. The scriptsassume local files prepared according to the manuscript Data Availability statement.

Reproducibility

Each active script writes a sessionInfo() or package-version file when run.The TCGA script also writes MD5 checksums for the public input files and records thegene-to-Ensembl mapping used.

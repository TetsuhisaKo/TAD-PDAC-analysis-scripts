# Data Dictionary

This document defines the input files, variables, generated intermediate files,
and public datasets required by the analysis scripts accompanying the manuscript:

**Clinically defined tumor-associated diabetes, intratumoral T-cell densities,
and survival in resected pancreatic ductal adenocarcinoma**

Only analyses that are included in the current manuscript or its current
Supplementary Tables/Figures are described below.

## Group labels and coding

In the clinical data, `TAD` is coded as `0 = non-TAD` and `1 = TAD`.
The four-category diabetes variable `dm_group` is coded as:

- `0 = No diabetes`
- `1 = Long-standing stable diabetes`
- `2 = New-onset diabetes`
- `3 = Worsening diabetes`

New-onset and worsening diabetes constitute TAD; no diabetes and long-standing
stable diabetes constitute non-TAD.

In the organoid analyses, the six TAD-derived lines are
`KYK019, KYK020, KYK067, KYK084, KYK090, KYK093`; the remaining 30 lines are
non-TAD.

## Data availability

Patient-level clinical data and patient-derived organoid data are not distributed
in this repository because of ethical and privacy restrictions.

Organoid RNA-sequencing and whole-genome sequencing data are available under
controlled access in the European Genome-phenome Archive under accession numbers
`EGAS00001007212` and `EGAS00001007211`, respectively.

TCGA-PAAD transcriptomic, clinical and survival data were obtained from the UCSC
Xena GDC Hub. PanCanAtlas ABSOLUTE tumor-purity estimates were obtained from the
Genomic Data Commons PanCanAtlas release
`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`
(GDC file UUID `4f277128-f793-4354-a13d-30cc7fe9f6b5`).

## Active scripts

| Script file | Analysis | Current manuscript mapping |
|---|---|---|
| `01_survival_analysis.R` | Kaplan–Meier, log-rank and Cox analyses for OS/CSS, diabetes phenotype and stage | Fig. 1; Supplementary Figs. S1–S2; Supplementary Tables S2–S4 |
| `02_deseq2_differential_expression.R` | Organoid RNA-seq differential expression and PCA | Fig. 3a,d,e; Supplementary Fig. S4; Supplementary Table S7 |
| `03_preranked_gsea.R` | Pre-ranked Hallmark and Reactome GSEA | Fig. 3b,c; Supplementary Table S6 |
| `04_secreted_factor_tpm.py` | TPM-based sensitivity analysis of selected genes | Supplementary Table S7 |
| `05_tcga_paad_analysis.R` | TCGA-PAAD GDF15, immune-gene, CYT, DDIT3, purity-adjusted and survival analyses | Fig. 4c–e; Supplementary Table S8; Supplementary Fig. S5 |
| `06_driver_gene_wgs_analysis.R` | Major PDAC driver-gene alteration comparison in organoids | Supplementary Fig. S3 |

---

## 1. `clinical_data.xlsx` — clinical cohort (n = 162)

Used by: **`01_survival_analysis.R`**

One row per patient.

| Column | Type | Coding / units | Definition |
|---|---|---|---|
| `os_time` | numeric | months | Overall survival time from surgery |
| `os_event` | integer 0/1 | 1 = death from any cause; 0 = censored | Overall-survival event indicator |
| `css_event` | integer 0/1 | 1 = PDAC death; 0 = censored | Cancer-specific-survival event indicator; non-PDAC deaths are censored |
| `TAD` | integer 0/1 | 0 = non-TAD; 1 = TAD | Tumor-associated diabetes status |
| `dm_group` | integer 0–3 | 0 = no DM; 1 = long-standing stable DM; 2 = new-onset DM; 3 = worsening DM | Four-category diabetes phenotype |
| `stage` | integer 0/1 | 0 = ≤IIA; 1 = ≥IIB | Pathologic stage group, UICC 8th edition |
| `margin` | integer 0/1 | 0 = R0; 1 = R1 | Resection-margin status |
| `NLR` | integer 0/1 | 0 = <2.1; 1 = ≥2.1 | Dichotomized neutrophil-to-lymphocyte ratio |
| `CA199` | integer 0/1 | 0 = ≤37; 1 = >37 U/mL | Dichotomized serum CA19-9 |
| `CA199_raw` | numeric | U/mL, >0 | Raw serum CA19-9 used for log2-continuous sensitivity analysis |
| `HbA1c` | numeric | % | Preoperative HbA1c used in sensitivity analyses |
| `adj_chemo` | integer 0/1 | 0 = No; 1 = Yes | Postoperative adjuvant chemotherapy |
| `neoadj` | integer 0/1 | 0 = No; 1 = Yes | Neoadjuvant chemotherapy |

### Current cohort checks

- total `n = 162`
- non-TAD `n = 98`; TAD `n = 64`
- no diabetes `n = 81`
- long-standing stable diabetes `n = 17`
- new-onset diabetes `n = 40`
- worsening diabetes `n = 24`
- stage ≥IIB `n = 119`
- neoadjuvant-excluded cohort `n = 144`
- OS deaths `n = 120`
- PDAC-specific deaths `n = 110`
- diabetes-only subset `n = 81`

Rows with `dm_group` 2 or 3 must have `TAD = 1`; rows with `dm_group` 0 or 1
must have `TAD = 0`.

The primary multivariable survival model includes TAD, pathologic stage,
resection-margin status, NLR and CA19-9. Sensitivity analyses additionally assess
continuous HbA1c, adjuvant chemotherapy, log2-continuous CA19-9, combined
log2-continuous CA19-9 plus adjuvant chemotherapy, and exclusion of patients who
received neoadjuvant chemotherapy.

The diabetes-only analysis compares TAD with long-standing stable diabetes among
the 81 patients with diabetes, with and without additional adjustment for
continuous HbA1c. The four-category model uses patients without diabetes as the
reference group. TAD-by-stage interaction models and stage-stratified adjusted
hazard ratios use the same clinical covariate framework.

`CA199_raw` is analyzed as `log2(CA199_raw)` without an offset.

---

## 2. `count_matrix.xlsx` — organoid RNA-seq raw counts (36 lines)

Used by: **`02_deseq2_differential_expression.R`**

| Column | Type | Definition |
|---|---|---|
| first column | string | Unique gene symbol used as the row name and for downstream MSigDB matching |
| one column per organoid line | integer | Raw read counts; column names begin with the organoid ID |

### Sample composition

- TAD: `n = 6`
- non-TAD: `n = 30`
- TAD lines: `KYK019, KYK020, KYK067, KYK084, KYK090, KYK093`
- MiSeq: `KYK015` and `KYK019`
- HiSeq 2500: all remaining lines

### DESeq2 specification

- primary design: `~ DM_group`
- reference group: non-TAD
- genes with row-summed raw counts `≤10` are excluded
- log2 fold changes are shrunken using `apeglm`
- differential expression is defined by BH-adjusted `P < 0.05`
- no fold-change threshold is imposed
- PCA uses variance-stabilized counts
- sequencing-platform sensitivity model: `~ platform + DM_group`

The first column must contain unique gene symbols. If the original quantification
uses Ensembl identifiers, these must be mapped to gene symbols before running the
script.

---

## 3. Primary DESeq2 output — generated intermediate file

Generated by: **`02_deseq2_differential_expression.R`**  
Used by: **`03_preranked_gsea.R`**

Recommended filename:

`results/organoid_RNAseq/DESeq2_primary_all_results.xlsx`

| Column | Type | Definition |
|---|---|---|
| `SYMBOL` | string | Gene symbol |
| `baseMean` | numeric | Mean normalized count |
| `raw_log2FoldChange` | numeric | Unshrunken DESeq2 log2 fold change, TAD versus non-TAD |
| `lfcSE_raw` | numeric | Standard error from the unshrunken DESeq2 result |
| `stat` | numeric | DESeq2 test statistic |
| `pvalue` | numeric | Nominal DESeq2 P value |
| `padj` | numeric | Benjamini–Hochberg-adjusted P value |
| `log2FoldChange` | numeric | `apeglm`-shrunken log2 fold change, TAD versus non-TAD |
| `lfcSE_shrunken` | numeric | Standard error associated with the shrunken effect estimate |

The nominal `pvalue` and `padj` are taken from the ordinary DESeq2 result; the
shrunken `log2FoldChange` is used for effect-size display and for the sign of the
GSEA ranking statistic.

The pre-ranked GSEA score is:

`sign(log2FoldChange) × -log10(pvalue)`

Supplementary Table S7 includes selected primary and sequencing-platform-adjusted
results and TPM-based sensitivity analyses.

---

## 4. `organoid_expression_TPM.xlsx` — organoid TPM matrix

Used by: **`04_secreted_factor_tpm.py`**

| Sheet | Column | Type | Definition |
|---|---|---|---|
| `TPM` | `gene_symbol` | string | Gene symbol |
| `TPM` | one `*_TPM` column per organoid line | numeric | TPM value; each column begins with the organoid ID and ends in `_TPM` |

The current TPM-based sensitivity analysis evaluates the selected genes displayed
in Fig. 3e. Group comparisons use `log2(TPM + 1)` values and a two-sided
Mann–Whitney U test.

The current manuscript specifically reports that GDF15 remained higher in
TAD-associated organoids in the TPM-based analysis (`P = 0.040`).

---

## 5. `organoid_driver_alterations.xlsx` — WGS-derived driver-gene calls

Used by: **`08_driver_gene_wgs_analysis.R`**

One row per organoid line.

| Column | Type | Coding | Definition |
|---|---|---|---|
| `SampleID` | string | e.g. `KYK019` | Organoid identifier |
| `KRAS_any` | integer 0/1 | 1 = altered | KRAS mutation |
| `TP53_any` | integer 0/1 | 1 = altered | TP53 mutation, heterozygous deletion, or copy-neutral LOH/UPD |
| `CDKN2A_any` | integer 0/1 | 1 = altered | CDKN2A mutation or homozygous deletion |
| `SMAD4_any` | integer 0/1 | 1 = altered | SMAD4 mutation, homozygous deletion, or heterozygous deletion |

These are gene-level alteration calls derived from the WGS analysis. The
repository script performs TAD-versus-non-TAD Fisher exact tests and does not
process raw WGS sequence data.

Current cohort composition is TAD `n = 6` and non-TAD `n = 30`.

---

## 6. TCGA-PAAD public inputs

Used by: **`07_tcga_paad_analysis.R`**

UCSC Xena GDC Hub data were accessed on **2026-04-28**. PanCanAtlas ABSOLUTE
purity data were accessed on **2026-07-31**.

| File | Required fields / definition |
|---|---|
| `TCGA-PAAD.star_fpkm-uq.tsv.gz` | Expression matrix; Ensembl gene IDs × TCGA samples |
| `TCGA-PAAD.clinical.tsv.gz` | Clinical annotations; requires `sample` and `primary_diagnosis.diagnoses` |
| `TCGA-PAAD.survival.tsv.gz` | Survival data; requires `sample`, `OS`, and `OS.time` |
| `gencode.v36.annotation.gtf.gene.probemap` | GENCODE v36 mapping; requires `id` and `gene` |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | PanCanAtlas ABSOLUTE purity; requires sample identifier, call status and purity |

### Expression scale

The manuscript analyzes TCGA expression on the log2-transformed FPKM-UQ scale.
The UCSC Xena matrix used for the final analysis already contains the transformed
values used in the analysis workflow. Therefore the values are used directly in
the repository script and **must not be log-transformed a second time**.

Ensembl version suffixes are removed after `_PAR_Y` records are excluded, and
gene identifiers are mapped using the GENCODE v36 probe map.

### Histology-restricted cohorts

The current expression cohort consists of 163 primary tumors with
PDAC-compatible histology:

- infiltrating duct carcinoma, NOS: `n = 143`
- adenocarcinoma, NOS: `n = 20`

Complete OS data are available for `n = 162`.

ABSOLUTE tumor-purity estimates are available for `n = 143` of the 163 tumors.

The histologically stringent sensitivity analysis is restricted to infiltrating
duct carcinoma, NOS:

- expression cohort: `n = 143`
- complete OS: `n = 142`
- ABSOLUTE purity available: `n = 123`

### GDF15 analyses

GDF15 is dichotomized at the median of the 163-tumor expression cohort. This
produces 81 GDF15-low and 82 GDF15-high tumors. After restriction to tumors with
complete OS data, the survival cohort contains 81 GDF15-low and 81 GDF15-high
tumors.

The six transcriptomic outcomes are:

`CD8A, GZMA, GZMB, PRF1, CYT, DDIT3`

CYT is calculated as the arithmetic mean of the log2-scale GZMA and PRF1
expression values.

The current analyses include:

1. two-sided Mann–Whitney U comparisons of the six outcomes between
   GDF15-low and GDF15-high tumors, with BH correction across these six tests;
2. continuous Spearman correlations between GDF15 and the same six outcomes;
3. partial Spearman correlations adjusted for ABSOLUTE tumor purity using
   rank-based residualization, with Fisher-z confidence intervals and BH
   correction across the six purity-adjusted correlations;
4. sensitivity analyses restricted to infiltrating duct carcinoma, NOS; and
5. Kaplan–Meier/log-rank analysis of OS according to median-dichotomized GDF15.

Current manuscript results include:

- GDF15 versus CD8A: Spearman `rho = -0.37`, `P < 0.001`
- GDF15 versus CYT: Spearman `rho = -0.40`, `P < 0.001`
- GDF15 versus DDIT3: Spearman `rho = 0.53`, `P < 0.001`
- GDF15 versus ABSOLUTE purity: Spearman `rho = 0.18`, `P = 0.031`
- purity-adjusted GDF15 versus CD8A: partial `rho = -0.33`
- purity-adjusted GDF15 versus CYT: partial `rho = -0.36`
- purity-adjusted GDF15 versus DDIT3: partial `rho = 0.51`
- GDF15-high versus GDF15-low OS: log-rank `P = 0.75`

---

## 7. GSEA output conventions

Used by: **`03_preranked_gsea.R`**

Pre-ranked GSEA uses:

- ranking statistic: `sign(shrunken log2FC) × -log10(nominal P)`
- MSigDB Hallmark and Reactome collections
- gene-set size: 15–500
- random seed: 42
- prespecified GSEA threshold: `q < 0.25`

Fig. 3b displays selected Hallmark signals and Fig. 3c displays selected Reactome
signals. Full pathway-level results are summarized in Supplementary Table S6.
Selected Hallmark signals that do not meet the prespecified FDR threshold are
interpreted as nominal signals rather than statistically significant enrichment.

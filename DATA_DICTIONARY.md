
# Data dictionary

This document describes the input files and variable coding required by the analysis scripts accompanying the TAD–PDAC manuscript. Patient-level clinical data and controlled-access organoid data are not included in this repository. All paths are relative to the repository root.

## 1. Surgical cohort: `data/clinical_data.xlsx`

Used by `01_survival_analysis.R` and `05_tcell_adjusted_models.R`. One row per patient (n = 162). The two scripts use the first worksheet.

### Required columns

| Column | Type / coding | Definition |
| --- | --- | --- |
| `os_time` | Numeric; months | Time from surgery to death or last follow-up; also used as the time variable for cancer-specific survival. |
| `os_event` | 0/1 | Death from any cause (`1`) or censored (`0`). |
| `css_event` | 0/1 | Death attributed to PDAC (`1`) or censored (`0`); deaths from other causes are censored. |
| `TAD` | 0/1 | Tumor-associated diabetes (`1`) versus non-TAD (`0`). |
| `dm_group` | 0/1/2/3 | `0`: no diabetes; `1`: long-standing stable diabetes; `2`: new-onset diabetes; `3`: worsening diabetes. |
| `stage` | 0/1 | Pathologic stage ≤IIA (`0`) or ≥IIB (`1`). |
| `margin` | 0/1 | R0 (`0`) or R1 (`1`) resection margin. |
| `NLR` | 0/1 | Neutrophil-to-lymphocyte ratio <2.1 (`0`) or ≥2.1 (`1`). |
| `NLR_raw` | Numeric; ≥0 | Continuous neutrophil-to-lymphocyte ratio. |
| `CA199` | 0/1 | Serum CA19-9 ≤37 U/mL (`0`) or >37 U/mL (`1`). |
| `CA199_raw` | Numeric; >0, U/mL | Continuous CA19-9, transformed as `log2(CA199_raw)` without an offset in sensitivity analyses. |
| `HbA1c` | Numeric; % | Preoperative glycated hemoglobin. |
| `adj_chemo` | 0/1 | Postoperative adjuvant chemotherapy: no (`0`), yes (`1`). |
| `neoadj` | 0/1 | Neoadjuvant chemotherapy: no (`0`), yes (`1`). |
| `age` | Numeric; years | Age at surgery. |
| `surgery_year` | Numeric | Calendar year of surgery. |
| `worsening_basis` | 1/2/3 or missing | For `dm_group = 3` only: `1` = HbA1c increase alone; `2` = treatment intensification alone; `3` = both. Missing for other groups. |
| `stable_longitudinal_data` | 0/1 or missing | For `dm_group = 1` only: `1` = longitudinal glycemic data available; `0` = not available. Missing for other groups. |
| `BMI` | Numeric; kg/m² | Body mass index. |
| `SMI` | Numeric; cm²/m² | Skeletal muscle index. |
| `albumin` | Numeric; g/dL | Serum albumin. |
| `ECOG_PS` | Integer | Eastern Cooperative Oncology Group performance status; categorized as 0 versus ≥1 for sensitivity analyses. |
| `CD4` | Numeric; >0, cells/mm² | Intratumoral CD4-positive T-cell density. |
| `CD8` | Numeric; >0, cells/mm² | Intratumoral CD8-positive T-cell density. |
| `FOXP3` | Numeric; >0, cells/mm² | Intratumoral FOXP3-positive cell density. |

`TAD = 1` corresponds to `dm_group` 2 or 3; `TAD = 0` corresponds to groups 0 or 1. The main Cox model uses the categorized `NLR` and `CA199` variables. T-cell density models use continuous `NLR_raw` and `log2(CA199_raw)`.

The survival script also uses `log2(CD8)` in its additional OS sensitivity model. Clinical data used by the main models must be complete; the two phenotype-specific classification fields may be missing outside their respective groups.

### Additional inputs for Table 2 univariable analyses

The following columns are accepted by `01_survival_analysis.R`. Each binary variable must contain complete 0/1 values, with `0` as its reference category. The accepted alternatives are listed in lookup order.

| Analysis term | Accepted input column(s) | Coding / units |
| --- | --- | --- |
| Age ≥65 years | `age` | Derived as age ≥65 versus <65 years. |
| Male sex | `male`, `male_sex`; alternatively `sex`, `Sex` | Binary: `0` = female, `1` = male; text alternative must contain `Female`/`Male`. Numeric `sex` coding is not inferred. |
| Tumor size >20 mm | `tumor_gt20`, `tumor_size_gt20`; alternatively `tumor_size_mm`, `tumor_size_mm_raw` | Binary: `0` = ≤20 mm, `1` = >20 mm; or size in mm. |
| KRAS mutation | `KRAS_mutation`, `KRAS_mut`, `KRAS_any` | `1` = mutation identified in the surgical cohort. |
| TP53 alteration | `TP53_alteration`, `TP53_alt`, `TP53_any` | `1` = altered according to the surgical-cohort assessment. |
| SMAD4 alteration | `SMAD4_alteration`, `SMAD4_alt`, `SMAD4_any` | `1` = altered according to the surgical-cohort assessment. |
| CDKN2A alteration | `CDKN2A_alteration`, `CDKN2A_alt`, `CDKN2A_any` | `1` = altered according to the surgical-cohort assessment. |
| CEA >5 ng/mL | `CEA_gt5`, `CEA_over5`; alternatively `CEA_raw`, `CEA_ng_mL` | Binary: `0` = ≤5 ng/mL, `1` = >5 ng/mL; or continuous CEA in ng/mL. |

The remaining Table 2 terms (`TAD`, `stage`, `margin`, `NLR`, `CA199`) use the required columns above, giving 13 univariable terms in total. The multivariable OS model uses these five prespecified terms only. If a required input is unavailable or invalid, the corresponding univariable estimate is reported as missing.

### Cohort validation

| Check | Expected value |
| --- | --- |
| Total surgical patients | 162 |
| Non-TAD / TAD | 98 / 64 |
| No / stable / new-onset / worsening diabetes | 81 / 17 / 40 / 24 |
| Stage ≤IIA / ≥IIB | 43 / 119 |
| OS / PDAC-specific deaths | 120 / 110 |
| Patients excluding neoadjuvant therapy | 144 |
| Worsening classification: HbA1c-only / treatment-only / both | 9 / 11 / 4 |
| Stable diabetes: longitudinal data available / unavailable | 8 / 9 |

## 2. Organoid RNA-seq counts: `data/count_matrix.xlsx`

Used by `02_deseq2_differential_expression.R`.

| Field | Format |
| --- | --- |
| First column | Unique gene symbols; one gene per row. |
| Remaining 36 columns | One column per organoid; nonnegative, integer raw read counts. |
| Column names | Begin with each organoid's sample ID; the script extracts the part preceding the first underscore. |

The cohort consists of 6 TAD-derived and 30 non-TAD-derived organoids. Two were sequenced on MiSeq and the remainder on HiSeq 2500. The current script determines group and sequencing platform from locally configured sample identifiers. Genes with summed raw counts ≤10 are excluded. The primary design is `~ DM_group`; the platform-adjusted design is `~ platform + DM_group`.

Log2 fold changes are shrunk using `apeglm`; Wald-test P values and BH-adjusted P values are taken from the unshrunken DESeq2 results. PCA uses variance-stabilized expression and the 500 most variable genes.

### Generated intermediate files

`02` writes `results/organoid_RNAseq/DESeq2_primary_all_results.xlsx`, used by `03` and `04`.

| Column | Definition |
| --- | --- |
| `SYMBOL` | Gene symbol. |
| `baseMean` | Mean normalized count. |
| `raw_log2FoldChange` | Unshrunken TAD-versus-non-TAD log2 fold change. |
| `lfcSE_raw` | Unshrunken log2 fold-change standard error. |
| `stat` | Wald statistic. |
| `pvalue` | Nominal Wald-test P value. |
| `padj` | BH-adjusted P value. |
| `log2FoldChange` | `apeglm`-shrunken log2 fold change. |
| `lfcSE_shrunken` | Shrunk estimate standard error. |

## 3. Gene-set enrichment analysis

`03_preranked_gsea.R` ranks genes from the DESeq2 output as:

`sign(log2FoldChange) × −log10(pvalue)`

The script uses MSigDB **2025.1.Hs** Hallmark and Reactome gene sets, 15–500 genes per set, with BH correction within each collection.

| Output | Reporting threshold |
| --- | --- |
| Supplementary Table S8A: Hallmark | BH-adjusted P <0.25 |
| Supplementary Table S8B / Fig. 3b: Reactome | BH-adjusted P <0.05 |
| Additional exploratory Reactome results | BH-adjusted P <0.25 |

Stage IV-adjusted analyses in Supplementary Table S8C require individual-level clinical annotations that are not distributed and are not implemented in the public scripts.

## 4. Organoid TPM matrix: `data/organoid_expression_TPM.xlsx`

Used by `04_secreted_factor_tpm.py`.

| Field | Format |
| --- | --- |
| Worksheet | `TPM` |
| `gene_symbol` | Gene identifier column. |
| Remaining expression columns | 36 organoid columns with names ending in `_TPM`; the sample ID precedes the first underscore. |
| Values | Nonnegative, finite TPM values. |

For the 12 genes shown in Fig. 3d, duplicate gene-symbol rows are averaged by sample. The script reports medians and quartiles on the TPM scale and two-sided, asymptotic Mann–Whitney U P values (continuity and tie corrections; no multiplicity correction). Primary DESeq2 estimates are merged from script `02` to produce Supplementary Table S9B.

## 5. Organoid driver-gene calls: `data/organoid_driver_alterations.xlsx`

Used by `07_driver_gene_wgs_analysis.R` to produce statistics for Supplementary Fig. S4. One row per organoid (36 rows).

| Column | Format / alteration definition |
| --- | --- |
| `SampleID` | Organoid sample identifier, matching the RNA-seq cohort. |
| `KRAS_any` | 0/1: KRAS mutation. |
| `TP53_any` | 0/1: TP53 mutation, heterozygous deletion, or copy-neutral loss of heterozygosity/uniparental disomy. |
| `CDKN2A_any` | 0/1: CDKN2A mutation or homozygous deletion. |
| `SMAD4_any` | 0/1: SMAD4 mutation, homozygous deletion, or heterozygous deletion. |

Calls are pre-derived from WGS. Script `07` performs two-sided Fisher exact tests without multiple-testing adjustment; it does not process reads or call variants.

## 6. TCGA-PAAD public inputs: `data/TCGA/`

Used by `06_tcga_paad_analysis.R`. The TCGA-PAAD files are obtained from public data resources.

| File | Required fields / content |
| --- | --- |
| `TCGA-PAAD.star_fpkm-uq.tsv.gz` | Gene-by-sample STAR FPKM-UQ expression matrix; primary tumors have sample IDs ending in `-01A`. |
| `TCGA-PAAD.clinical.tsv.gz` | Clinical annotations; `sample`, `primary_diagnosis.diagnoses`. |
| `TCGA-PAAD.survival.tsv.gz` | Survival annotations; `sample`, `OS`, `OS.time`. |
| `gencode.v36.annotation.gtf.gene.probemap` | Gene map; `id`, `gene`. |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | PanCanAtlas ABSOLUTE purity; `array`, `call status`, `purity`. |

UCSC Xena GDC Hub expression data (accessed 28 April 2026) are used on their supplied **log2(FPKM-UQ + 1)** scale; do not transform them again. Gene symbols are mapped with GENCODE v36. The expression cohort has 163 PDAC-compatible primary tumors (143 infiltrating duct carcinoma, NOS; 20 adenocarcinoma, NOS), 162 with OS data and 143 with ABSOLUTE purity estimates. The ductal-histology-only cohort contains 143 tumors (123 with purity estimates).

GDF15 is split at the median of the 163 tumors (82 high, 81 low). The transcriptomic outcomes are `CD8A`, `GZMA`, `GZMB`, `PRF1`, `CYT`, and `DDIT3`. `CYT` is the arithmetic mean of log2-scale `GZMA` and `PRF1` expression. TAD phenotype is not defined in TCGA-PAAD because longitudinal diabetes-course information is unavailable.

## Data access

Clinical and organoid patient-level inputs are not publicly distributed. Organoid RNA-seq and WGS sequence data are available under controlled access through the European Genome-phenome Archive, studies `EGAS00001007212` and `EGAS00001007211`, respectively. All donor-level input files and generated files containing donor-linked metadata should remain outside the public repository.


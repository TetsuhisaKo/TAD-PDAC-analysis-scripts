
#!/usr/bin/env python3
"""
04_secreted_factor_tpm.py

Outputs (results/organoid_TPM/):
  Suppl. Table S9B  TPM-based sensitivity analysis for the
                    Fig. 3d genes (with DESeq2 estimates from
                    02_deseq2_differential_expression.R)
  Fig. 3d shows DESeq2 log2FC (script 02).

- 36 organoids: TAD n=6, non-TAD n=30
- Two-sided Mann-Whitney U test, normal approximation with
  continuity and tie correction. The test is rank-based, so
  P values are identical on the TPM and log2(TPM + 1) scales.
- Median/IQR on the original TPM scale (linear interpolation)
- No multiplicity adjustment for TPM sensitivity P values
- DESeq2 columns: apeglm-shrunken log2FC and BH-adjusted P
  (padj) from the primary model
"""

from pathlib import Path
import sys
import warnings
import numpy as np
import pandas as pd
import scipy
from scipy import stats

INPUT_PATH = Path("data/organoid_expression_TPM.xlsx")
SHEET = "TPM"
DESEQ_PATH = Path(
    "results/organoid_RNAseq/DESeq2_primary_all_results.xlsx"
)
META_PATH = Path(
    "results/organoid_RNAseq/organoid_sample_metadata.csv"
)
OUT_DIR = Path("results/organoid_TPM")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENE_COL = "gene_symbol"

TAD_CASES = {
    "KYK019", "KYK020", "KYK067",
    "KYK084", "KYK090", "KYK093",
}

# Fig. 3d genes (same as fig3d_genes in 02)
GENES = [
    "GDF15", "PRSS2", "CPZ", "NXPH4",
    "LYZ", "P3H3", "HHIPL1", "SCG2",
    "IL33", "COPA", "BMP4", "FGF19",
]

if not INPUT_PATH.exists():
    raise FileNotFoundError(f"Missing input: {INPUT_PATH}")

df = pd.read_excel(INPUT_PATH, sheet_name=SHEET)

if GENE_COL not in df.columns:
    raise ValueError(f"Required column {GENE_COL} not found.")

tpm_cols = [
    c for c in df.columns
    if str(c).endswith("_TPM")
]

if len(tpm_cols) != 36:
    raise ValueError(
        f"Expected 36 TPM columns, got {len(tpm_cols)}."
    )

sample_ids = [str(c).split("_")[0] for c in tpm_cols]

if len(set(sample_ids)) != len(sample_ids):
    raise ValueError("Duplicate organoid IDs among TPM columns.")
if not TAD_CASES.issubset(sample_ids):
    missing = sorted(TAD_CASES - set(sample_ids))
    raise ValueError(f"TAD organoids not found: {missing}")

groups = np.array([
    "TAD" if s in TAD_CASES else "non-TAD"
    for s in sample_ids
])

if np.sum(groups == "TAD") != 6:
    raise ValueError("TAD organoid count must be 6.")
if np.sum(groups == "non-TAD") != 30:
    raise ValueError("Non-TAD organoid count must be 30.")

# Cross-check samples and groups against the DESeq2 metadata
if META_PATH.exists():
    meta = pd.read_csv(META_PATH)
    meta_groups = dict(zip(meta["case"], meta["DM_group"]))
    if set(meta_groups) != set(sample_ids):
        raise ValueError(
            "Organoids in the TPM file differ from the DESeq2 "
            "count matrix."
        )
    mismatch = [
        s for s, g in zip(sample_ids, groups)
        if meta_groups[s] != g
    ]
    if mismatch:
        raise ValueError(f"Group assignment differs: {mismatch}")
else:
    warnings.warn(
        f"{META_PATH} not found; sample cross-check skipped."
    )


def q1(x):
    return float(np.quantile(x, 0.25))


def q3(x):
    return float(np.quantile(x, 0.75))


rows = []

for gene in GENES:
    sub = df.loc[df[GENE_COL] == gene, tpm_cols]

    if sub.empty:
        raise ValueError(f"Gene not found: {gene}")

    # Duplicate gene-symbol rows are collapsed by sample mean.
    if len(sub) > 1:
        warnings.warn(
            f"{gene}: {len(sub)} rows collapsed by sample mean."
        )
    raw_tpm = sub.astype(float).mean(axis=0).to_numpy()

    if np.any(~np.isfinite(raw_tpm)):
        raise ValueError(f"Non-finite TPM for {gene}")
    if np.any(raw_tpm < 0):
        raise ValueError(f"Negative TPM for {gene}")

    non_raw = raw_tpm[groups == "non-TAD"]
    tad_raw = raw_tpm[groups == "TAD"]

    test = stats.mannwhitneyu(
        tad_raw,
        non_raw,
        alternative="two-sided",
        method="asymptotic",
        use_continuity=True,
    )

    rows.append({
        "Gene": gene,
        "Rows_collapsed": len(sub),
        "Non_TAD_n": len(non_raw),
        "TAD_n": len(tad_raw),
        "Non_TAD_TPM_median": float(np.median(non_raw)),
        "Non_TAD_TPM_Q1": q1(non_raw),
        "Non_TAD_TPM_Q3": q3(non_raw),
        "TAD_TPM_median": float(np.median(tad_raw)),
        "TAD_TPM_Q1": q1(tad_raw),
        "TAD_TPM_Q3": q3(tad_raw),
        "TPM_based_P": float(test.pvalue),
    })

out = pd.DataFrame(rows)

# DESeq2 estimates from the primary model
if DESEQ_PATH.exists():
    deseq = pd.read_excel(DESEQ_PATH)
    required = {"SYMBOL", "log2FoldChange", "padj"}
    if not required.issubset(deseq.columns):
        raise ValueError(
            f"{DESEQ_PATH} lacks columns: "
            f"{sorted(required - set(deseq.columns))}"
        )
    deseq = (
        deseq.loc[
            deseq["SYMBOL"].isin(GENES),
            ["SYMBOL", "log2FoldChange", "padj"],
        ]
        .rename(columns={
            "SYMBOL": "Gene",
            "log2FoldChange": "DESeq2_log2FC",
            "padj": "DESeq2_padj",
        })
    )
    out = out.merge(deseq, on="Gene", how="left")

    absent = out.loc[out["DESeq2_log2FC"].isna(), "Gene"].tolist()
    no_padj = out.loc[
        out["DESeq2_log2FC"].notna() & out["DESeq2_padj"].isna(),
        "Gene",
    ].tolist()
    if absent:
        warnings.warn(
            f"Not in DESeq2 results (count filter): {absent}"
        )
    if no_padj:
        warnings.warn(
            f"DESeq2 padj missing (independent filtering or "
            f"Cook's outlier): {no_padj}"
        )
else:
    warnings.warn(
        f"{DESEQ_PATH} not found; DESeq2 columns omitted. "
        "Run 02_deseq2_differential_expression.R first."
    )

preferred_cols = [
    "Gene",
    "DESeq2_log2FC",
    "DESeq2_padj",
    "Non_TAD_TPM_median",
    "Non_TAD_TPM_Q1",
    "Non_TAD_TPM_Q3",
    "TAD_TPM_median",
    "TAD_TPM_Q1",
    "TAD_TPM_Q3",
    "TPM_based_P",
    "Non_TAD_n",
    "TAD_n",
    "Rows_collapsed",
]

out = out[
    [col for col in preferred_cols if col in out.columns]
]

csv_path = OUT_DIR / "Supplementary_Table_S9B_TPM_sensitivity.csv"
xlsx_path = OUT_DIR / "Supplementary_Table_S9B_TPM_sensitivity.xlsx"

out.to_csv(csv_path, index=False)
out.to_excel(xlsx_path, index=False)

print(out.to_string(index=False))

version_path = OUT_DIR / "python_package_versions.txt"
version_path.write_text(
    f"python={sys.version.split()[0]}\n"
    f"pandas={pd.__version__}\n"
    f"numpy={np.__version__}\n"
    f"scipy={scipy.__version__}\n"
    "test=scipy.stats.mannwhitneyu(two-sided, asymptotic, "
    "continuity correction)\n",
    encoding="utf-8",
)

print("\nTPM sensitivity analyses completed.")

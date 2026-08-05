# ============================================================
# 04_secreted_factor_tpm.py
#
# TPM-based sensitivity analysis of selected genes displayed
# in Fig. 3e.
#
# Computes Supplementary Table S7B.
# Plotting code is intentionally not included.
#
# TAD n=6 vs non-TAD n=30.
# Two-sided Mann-Whitney U, asymptotic method.
#
# The test is performed on log2(TPM+1). Because log2(TPM+1)
# is strictly monotonic, ranks and ties are unchanged, so the
# Mann-Whitney P value is identical to that obtained from raw
# TPM. The transformed scale is retained for consistency with
# the expression-scale description.
#
# TPM-based P values are nominal/unadjusted.
# ============================================================

from pathlib import Path
import sys
import pandas as pd
import numpy as np
import scipy
from scipy import stats

INPUT_PATH = Path("data/organoid_expression_TPM.xlsx")
DESEQ_PATH = Path("results/organoid_RNAseq/DESeq2_primary_all_results.xlsx")
OUT_DIR = Path("results/organoid_TPM")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SHEET = "TPM"
GENE_COL = "gene_symbol"

if not INPUT_PATH.exists():
    raise FileNotFoundError(f"Missing input file: {INPUT_PATH}")

df_expr = pd.read_excel(INPUT_PATH, sheet_name=SHEET)

if GENE_COL not in df_expr.columns:
    raise ValueError(f"'{GENE_COL}' column not found.")

tpm_cols = [c for c in df_expr.columns if str(c).endswith("_TPM")]

if len(tpm_cols) != 36:
    raise ValueError(
        f"Expected 36 *_TPM sample columns, found {len(tpm_cols)}."
    )

TAD_CASES = {
    "KYK019", "KYK020", "KYK067",
    "KYK084", "KYK090", "KYK093",
}

sample_ids = [str(c).split("_")[0] for c in tpm_cols]
groups = np.array([
    "TAD" if s in TAD_CASES else "non-TAD"
    for s in sample_ids
])

if int(np.sum(groups == "TAD")) != 6:
    raise ValueError("Expected 6 TAD organoid columns.")
if int(np.sum(groups == "non-TAD")) != 30:
    raise ValueError("Expected 30 non-TAD organoid columns.")

# Current Fig. 3e selected genes.
GENES = [
    "GDF15", "PRSS2", "CPZ", "NXPH4", "LYZ", "P3H3",
    "HHIPL1", "SCG2", "IL33", "COPA", "BMP4", "FGF19",
]

def q1(x):
    return float(np.quantile(x, 0.25))

def q3(x):
    return float(np.quantile(x, 0.75))

rows = []

for gene in GENES:
    gene_rows = df_expr.loc[
        df_expr[GENE_COL].astype(str) == gene,
        tpm_cols
    ]

    if gene_rows.empty:
        raise ValueError(f"Gene not found in TPM matrix: {gene}")

    # If duplicate gene-symbol rows exist, collapse sample-wise by mean.
    values = gene_rows.astype(float).mean(axis=0).to_numpy()

    if np.any(~np.isfinite(values)):
        raise ValueError(f"Non-finite TPM value detected for {gene}")
    if np.any(values < 0):
        raise ValueError(f"Negative TPM value detected for {gene}")

    non_raw = values[groups == "non-TAD"]
    tad_raw = values[groups == "TAD"]

    non_log = np.log2(non_raw + 1.0)
    tad_log = np.log2(tad_raw + 1.0)

    test = stats.mannwhitneyu(
        tad_log,
        non_log,
        alternative="two-sided",
        method="asymptotic",
        use_continuity=True,
    )

    rows.append({
        "Gene": gene,
        "Non_TAD_n": int(len(non_raw)),
        "TAD_n": int(len(tad_raw)),
        "Non_TAD_TPM_median": float(np.median(non_raw)),
        "Non_TAD_TPM_Q1": q1(non_raw),
        "Non_TAD_TPM_Q3": q3(non_raw),
        "TAD_TPM_median": float(np.median(tad_raw)),
        "TAD_TPM_Q1": q1(tad_raw),
        "TAD_TPM_Q3": q3(tad_raw),
        "TPM_based_P": float(test.pvalue),
    })

out = pd.DataFrame(rows)

# Merge primary DESeq2 results when available.
if DESEQ_PATH.exists():
    deseq = pd.read_excel(DESEQ_PATH)
    required = {"SYMBOL", "log2FoldChange", "padj"}

    if required.issubset(deseq.columns):
        deseq = (
            deseq.loc[
                deseq["SYMBOL"].isin(GENES),
                ["SYMBOL", "log2FoldChange", "padj"],
            ]
            .rename(columns={
                "SYMBOL": "Gene",
                "log2FoldChange": "DESeq2_log2FC",
                "padj": "DESeq2_adjusted_P",
            })
        )
        out = out.merge(deseq, on="Gene", how="left")

order_map = {g: i for i, g in enumerate(GENES)}
out["_order"] = out["Gene"].map(order_map)
out = out.sort_values("_order").drop(columns="_order")

preferred_cols = [
    "Gene",
    "DESeq2_log2FC",
    "DESeq2_adjusted_P",
    "Non_TAD_TPM_median",
    "Non_TAD_TPM_Q1",
    "Non_TAD_TPM_Q3",
    "TAD_TPM_median",
    "TAD_TPM_Q1",
    "TAD_TPM_Q3",
    "TPM_based_P",
    "Non_TAD_n",
    "TAD_n",
]

out = out[[c for c in preferred_cols if c in out.columns]]

out.to_csv(
    OUT_DIR / "Supplementary_Table_S7B_TPM_sensitivity.csv",
    index=False
)

out.to_excel(
    OUT_DIR / "Supplementary_Table_S7B_TPM_sensitivity.xlsx",
    index=False
)

print(out.to_string(index=False))

for gene, expected in [("GDF15", "~0.040"), ("IL33", "~0.155")]:
    row = out.loc[out["Gene"] == gene]
    if not row.empty:
        print(
            f"{gene} manuscript check: "
            f"P={float(row['TPM_based_P'].iloc[0]):.6g} "
            f"(expected {expected})"
        )

(OUT_DIR / "python_package_versions.txt").write_text(
    f"python={sys.version.split()[0]}\n"
    f"pandas={pd.__version__}\n"
    f"numpy={np.__version__}\n"
    f"scipy={scipy.__version__}\n",
    encoding="utf-8",
)

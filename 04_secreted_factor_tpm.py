# ============================================================
# Script 3: Secreted-factor comparison, TAD vs non-TAD (TPM)
#   Two-sided Mann-Whitney U on log2(TPM + 1)
#
# Maps to: GDF15/IL33 TPM-based sensitivity analysis shown in Suppl. Fig. S5.
#          Extend the `genes` list to reproduce the full secreted-factor panel.
# Group labels: "DM" = tumor-associated diabetes (TAD); "Normal" = non-TAD.
#
# This TPM-based comparison is a normalization-independent sensitivity
# analysis complementing the DESeq2 (raw-count) results. It reproduces the
# manuscript values GDF15 p = 0.040 and IL33 p = 0.15.
#
# Reproducibility note: scipy.stats.mannwhitneyu selects its default method
# from sample size / ties, and this default has changed across SciPy versions
# (newer versions may select the exact test). The manuscript values were
# computed with the normal approximation, so method="asymptotic" is set
# explicitly here to reproduce them regardless of SciPy version.
#
# Input data:
#   data/organoid_expression_TPM.xlsx
#     - This file is NOT included in the repository: the underlying
#       patient-derived organoid transcriptomic data are under controlled
#       access and subject to ethical restrictions. Availability follows the
#       Data Availability statement of the manuscript.
#     - Users with appropriate access should prepare an input file with a
#       sheet named "TPM", a "gene_symbol" column, and one *_TPM column per
#       organoid line (sample column names ending in "_TPM"). The six
#       TAD ("DM") lines are listed in DM_CASES below.
#
# Usage:
#   python 04_secreted_factor_tpm.py
#
# Tested with: Python 3.10.12, pandas 2.2.3, numpy 2.2.5, scipy 1.15.3
# ============================================================
import sys
import pandas as pd
import numpy as np
import scipy
from scipy import stats

# ---- Input (relative path within the repository; data not distributed here) ----
INPUT_PATH = "data/organoid_expression_TPM.xlsx"
SHEET = "TPM"

df_expr = pd.read_excel(INPUT_PATH, sheet_name=SHEET)

# ---- Identify gene-symbol column and TPM sample columns ----
GENE_COL = "gene_symbol"
tpm_cols = [c for c in df_expr.columns if str(c).endswith("_TPM")]  # 36 organoid lines
assert GENE_COL in df_expr.columns, f"'{GENE_COL}' column not found."
assert len(tpm_cols) > 0, "No *_TPM sample columns found."

# ---- Group assignment (leading sample ID before the first underscore) ----
# TAD-derived organoid lines (n = 6); all remaining lines are non-TAD (n = 30).
DM_CASES = {"KYK019", "KYK020", "KYK067", "KYK084", "KYK090", "KYK093"}
meta = pd.DataFrame({
    "Sample": tpm_cols,
    "Group":  ["DM" if str(c).split("_")[0] in DM_CASES else "Normal"
               for c in tpm_cols],
})
print(f"Samples: {len(tpm_cols)} total "
      f"(DM={int((meta.Group=='DM').sum())}, Normal={int((meta.Group=='Normal').sum())})")

# ---- Candidate secreted factors (extend to the full panel for Suppl. Fig. S5) ----
genes = ["GDF15", "IL33"]

for gene in genes:
    # Collapse duplicate gene symbols by mean.
    row = (df_expr[df_expr[GENE_COL] == gene]
           .groupby(GENE_COL)[tpm_cols].mean())
    if row.empty:
        print(f"{gene}: not found in {GENE_COL}")
        continue
    expr = pd.DataFrame({
        "Sample":  tpm_cols,
        "log2TPM": np.log2(row.values.flatten() + 1),
    }).merge(meta, on="Sample")

    dm_g  = expr.loc[expr["Group"] == "DM",     "log2TPM"]
    nor_g = expr.loc[expr["Group"] == "Normal", "log2TPM"]
    # method="asymptotic" reproduces the manuscript p-values across SciPy versions.
    mw    = stats.mannwhitneyu(dm_g, nor_g, alternative="two-sided",
                               method="asymptotic")
    lfc   = dm_g.mean() - nor_g.mean()

    print(f"{gene}: log2FC={lfc:.3f}, "
          f"DM median={dm_g.median():.3f}, "
          f"Normal median={nor_g.median():.3f}, "
          f"p={mw.pvalue:.4f}")

# ---- Record package versions for reproducibility ----
print(f"\n# python={sys.version.split()[0]}, pandas={pd.__version__}, "
      f"numpy={np.__version__}, scipy={scipy.__version__}")

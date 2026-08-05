"""
04_secreted_factor_tpm.py

TPM-based sensitivity analysis of selected genes shown in Fig. 3e.

Current manuscript / Supplementary Table S7B:
- 36 patient-derived PDAC organoids
- TAD n=6; non-TAD n=30
- Two-sided Mann-Whitney U test on log2(TPM + 1)
- Asymptotic approximation is specified explicitly
- Descriptive median and IQR are reported on the ORIGINAL TPM scale
- No multiple-testing correction is applied to the TPM-based sensitivity P values

The 12 current Fig. 3e genes are:
GDF15, PRSS2, CPZ, NXPH4, LYZ, P3H3, HHIPL1, SCG2,
IL33, COPA, BMP4, FGF19.

Optional:
If results/organoid_RNAseq/DESeq2_primary_all_results.xlsx exists,
DESeq2 shrunken log2FC and adjusted P are merged into the output.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats

INPUT_PATH = Path("data/organoid_expression_TPM.xlsx")
SHEET = "TPM"
DESEQ_PATH = Path("results/organoid_RNAseq/DESeq2_primary_all_results.xlsx")
OUT_DIR = Path("results/organoid_TPM")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENE_COL = "gene_symbol"

TAD_CASES = {
    "KYK019", "KYK020", "KYK067",
    "KYK084", "KYK090", "KYK093",
}

GENES = [
    "GDF15", "PRSS2", "CPZ", "NXPH4",
    "LYZ", "P3H3", "HHIPL1", "SCG2",
    "IL33", "COPA", "BMP4", "FGF19",
]

if not INPUT_PATH.exists():
    raise FileNotFoundError(
        f"{INPUT_PATH} not found. "
        "This controlled-access input is not distributed with the repository."
    )

df = pd.read_excel(INPUT_PATH, sheet_name=SHEET)

if GENE_COL not in df.columns:
    raise ValueError(f"Required column '{GENE_COL}' not found.")

tpm_cols = [c for c in df.columns if str(c).endswith("_TPM")]
if len(tpm_cols) != 36:
    raise ValueError(f"Expected 36 *_TPM columns, found {len(tpm_cols)}.")

sample_cases = [str(c).split("_")[0] for c in tpm_cols]
groups = np.array(["TAD" if c in TAD_CASES else "non-TAD" for c in sample_cases])

if int(np.sum(groups == "TAD")) != 6 or int(np.sum(groups == "non-TAD")) != 30:
    raise ValueError(
        "Group counts do not match the current cohort "
        f"(TAD={np.sum(groups=='TAD')}, non-TAD={np.sum(groups=='non-TAD')})."
    )

def q1(x):
    return float(np.quantile(np.asarray(x, dtype=float), 0.25))

def q3(x):
    return float(np.quantile(np.asarray(x, dtype=float), 0.75))

rows = []

for gene in GENES:
    sub = df.loc[df[GENE_COL] == gene, tpm_cols]

    if sub.empty:
        raise ValueError(f"{gene} not found in {GENE_COL}.")

    # If duplicate gene-symbol rows exist, collapse by arithmetic mean per sample.
    raw_tpm = sub.astype(float).mean(axis=0).to_numpy()

    non_raw = raw_tpm[groups == "non-TAD"]
    tad_raw = raw_tpm[groups == "TAD"]

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

# Optional merge with the current primary DESeq2 results.
if DESEQ_PATH.exists():
    deseq = pd.read_excel(DESEQ_PATH)
    required = {"SYMBOL", "log2FoldChange", "padj"}
    if required.issubset(deseq.columns):
        deseq = (
            deseq.loc[deseq["SYMBOL"].isin(GENES),
                      ["SYMBOL", "log2FoldChange", "padj"]]
            .rename(columns={
                "SYMBOL": "Gene",
                "log2FoldChange": "DESeq2_log2FC",
                "padj": "DESeq2_adjusted_P",
            })
        )
        out = out.merge(deseq, on="Gene", how="left")

# Arrange to mirror Supplementary Table S7B.
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

csv_path = OUT_DIR / "Supplementary_Table_S7B_TPM_sensitivity.csv"
xlsx_path = OUT_DIR / "Supplementary_Table_S7B_TPM_sensitivity.xlsx"

out.to_csv(csv_path, index=False)
out.to_excel(xlsx_path, index=False)

print(out.to_string(index=False))

# Current manuscript checks:
# GDF15: non-TAD median ~202.645; TAD median ~574.450; P ~0.040
# IL33:  non-TAD median ~0.425;   TAD median ~1.655;   P ~0.155

version_path = OUT_DIR / "python_package_versions.txt"
version_path.write_text(
    f"python={sys.version.split()[0]}\n"
    f"pandas={pd.__version__}\n"
    f"numpy={np.__version__}\n"
    f"scipy={scipy.__version__}\n",
    encoding="utf-8",
)sion__}")

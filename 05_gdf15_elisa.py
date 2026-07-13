# ============================================================
# ELISA-1: GDF15 ELISA quantification by 4PL standard curve
#
# Purpose:
#   Fit a four-parameter logistic (4PL) standard curve to GDF15 ELISA
#   absorbance data, interpolate GDF15 concentrations in conditioned
#   media, and normalize secreted GDF15 to organoid DNA content
#   as pg GDF15 per ug DNA.
#
# Maps to:
#   ELISA quantification used for Supplementary Fig. S6.
#
# Input data:
#   data/GDF15_ELISA_raw.xlsx
#     - This file is NOT included in the repository because the underlying
#       patient-derived organoid data are subject to ethical restrictions.
#     - Required sheets:
#         standard_curve: concentration_pg_ml, A450, A570
#         samples: Sample, A450, A570, dilution_factor,
#                  conditioned_medium_volume_ml
#         dna: Sample, DNA_ug
#
# Notes:
#   - OD is corrected as A450 - A570.
#   - The mean zero-standard OD is used as blank.
#   - Zero-standard wells are used for blank subtraction but excluded
#     from the 4PL fitting.
#   - conditioned_medium_volume_ml should represent the total recovered
#     conditioned-medium volume used to calculate total secreted GDF15.
#
# Usage:
#   python3 05_gdf15_elisa.py
#
# Tested with: Python 3.10.12, pandas 2.2.3, numpy 2.2.5, scipy 1.15.3
# ============================================================

import sys
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import curve_fit
from pathlib import Path

# ---- Paths ----
INPUT_PATH = Path("data/GDF15_ELISA_raw.xlsx")
OUTPUT_DIR = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)

OUTPUT_PATH = OUTPUT_DIR / "GDF15_ELISA_normalized.xlsx"

if not INPUT_PATH.exists():
    raise FileNotFoundError(
        f"Input file not found: {INPUT_PATH}. "
        "Prepare the input file according to the required sheet structure."
    )

# ---- Load input sheets ----
std = pd.read_excel(INPUT_PATH, sheet_name="standard_curve")
sam = pd.read_excel(INPUT_PATH, sheet_name="samples")
dna = pd.read_excel(INPUT_PATH, sheet_name="dna")

# ---- Basic column checks ----
required_std = {"concentration_pg_ml", "A450", "A570"}
required_sam = {"Sample", "A450", "A570", "dilution_factor",
                "conditioned_medium_volume_ml"}
required_dna = {"Sample", "DNA_ug"}

missing_std = required_std - set(std.columns)
missing_sam = required_sam - set(sam.columns)
missing_dna = required_dna - set(dna.columns)

if missing_std:
    raise ValueError(f"Missing columns in standard_curve sheet: {missing_std}")
if missing_sam:
    raise ValueError(f"Missing columns in samples sheet: {missing_sam}")
if missing_dna:
    raise ValueError(f"Missing columns in dna sheet: {missing_dna}")

# ---- OD correction and blank subtraction ----
for dat in (std, sam):
    dat["OD"] = dat["A450"].astype(float) - dat["A570"].astype(float)

blank_od = std.loc[std["concentration_pg_ml"] == 0, "OD"].mean()
if np.isnan(blank_od):
    blank_od = 0.0

std["OD_blank_corrected"] = std["OD"] - blank_od
sam["OD_blank_corrected"] = sam["OD"] - blank_od

# ---- Fit increasing 4PL curve ----
# 4PL model:
#   y = bottom + amplitude / (1 + (ec50 / x)^hill)
# where amplitude > 0, ec50 > 0, and hill > 0.
# This parameterization enforces an increasing curve for sandwich ELISA data.

def four_pl(x, bottom, amplitude, ec50, hill):
    x = np.asarray(x, dtype=float)
    return bottom + amplitude / (1.0 + (ec50 / x) ** hill)

fit_dat = std.loc[std["concentration_pg_ml"] > 0].copy()

x = fit_dat["concentration_pg_ml"].astype(float).to_numpy()
y = fit_dat["OD_blank_corrected"].astype(float).to_numpy()

p0 = [
    np.min(y),                  # bottom
    np.max(y) - np.min(y),       # amplitude
    np.median(x),                # ec50
    1.0                          # hill
]

bounds = (
    [-np.inf, 1e-12, 1e-12, 1e-12],
    [ np.inf, np.inf, np.inf, np.inf]
)

pars, cov = curve_fit(
    four_pl,
    x,
    y,
    p0=p0,
    bounds=bounds,
    maxfev=100000
)

bottom, amplitude, ec50, hill = pars
top = bottom + amplitude

print("4PL parameters:")
print({
    "bottom": bottom,
    "top": top,
    "amplitude": amplitude,
    "ec50": ec50,
    "hill": hill
})

# ---- Invert 4PL to interpolate concentration ----
def inverse_four_pl(y_obs, bottom, amplitude, ec50, hill):
    y_obs = np.asarray(y_obs, dtype=float)

    # From:
    #   y = bottom + amplitude / (1 + (ec50 / x)^hill)
    # Therefore:
    #   x = ec50 / ((amplitude / (y - bottom) - 1)^(1/hill))
    ratio = amplitude / (y_obs - bottom) - 1.0

    with np.errstate(divide="ignore", invalid="ignore"):
        conc = ec50 / (ratio ** (1.0 / hill))

    # Invalid values occur when OD is outside the valid fitted range.
    conc[~np.isfinite(conc)] = np.nan
    conc[ratio <= 0] = np.nan

    return conc

# ---- Summarize sample wells ----
sample_summary = (
    sam.groupby("Sample", as_index=False)
    .agg(
        OD_mean=("OD_blank_corrected", "mean"),
        OD_sd=("OD_blank_corrected", "std"),
        dilution_factor=("dilution_factor", "first"),
        conditioned_medium_volume_ml=("conditioned_medium_volume_ml", "first")
    )
)

# ---- Flag values outside the observed fitted standard-curve OD range ----
od_min, od_max = np.min(y), np.max(y)

sample_summary["outside_standard_range"] = (
    (sample_summary["OD_mean"] < od_min) |
    (sample_summary["OD_mean"] > od_max)
)

# ---- Interpolate concentrations ----
sample_summary["GDF15_pg_ml_raw"] = inverse_four_pl(
    sample_summary["OD_mean"].to_numpy(),
    bottom,
    amplitude,
    ec50,
    hill
)

sample_summary["GDF15_pg_ml"] = (
    sample_summary["GDF15_pg_ml_raw"] *
    sample_summary["dilution_factor"].fillna(1.0)
)

sample_summary["invalid_interpolation"] = sample_summary["GDF15_pg_ml_raw"].isna()

# ---- Normalize to DNA content ----
out = sample_summary.merge(
    dna[["Sample", "DNA_ug"]],
    on="Sample",
    how="left"
)

out["GDF15_pg_per_ug_DNA"] = (
    out["GDF15_pg_ml"] *
    out["conditioned_medium_volume_ml"].fillna(1.0) /
    out["DNA_ug"]
)

# ---- Save output ----
out.to_excel(OUTPUT_PATH, index=False)

print("\nNormalized GDF15 output:")
print(
    out[[
        "Sample",
        "OD_mean",
        "GDF15_pg_ml",
        "DNA_ug",
        "GDF15_pg_per_ug_DNA",
        "outside_standard_range",
        "invalid_interpolation"
    ]]
)

print(f"\nSaved: {OUTPUT_PATH}")

# ---- Record package versions for reproducibility ----
print(
    f"\n# python={sys.version.split()[0]}, "
    f"pandas={pd.__version__}, "
    f"numpy={np.__version__}, "
    f"scipy={scipy.__version__}"
)

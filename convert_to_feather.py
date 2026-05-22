"""Convert XP and CRISPR RDS files to feather for dkg ingestion."""

import pyreadr
import pandas as pd

DATA_DIR = "C:/GitHub/DepMap/data/26Q1"
META_COLS = {"X", "SequencingID", "ModelConditionID", "ModelID",
             "IsDefaultEntryForMC", "IsDefaultEntryForModel"}

# --- CRISPR ---
print("Reading CRISPR_26Q1.rds...")
crispr = pyreadr.read_r(f"{DATA_DIR}/CRISPR_26Q1.rds")[None]
# "X" column is the ModelID row identifier
crispr = crispr.rename(columns={"X": "row_id"})
print(f"  CRISPR: {crispr.shape[0]} rows x {crispr.shape[1]-1} gene cols")
print(f"  row_id sample: {crispr['row_id'].head(3).tolist()}")
crispr.to_feather(f"{DATA_DIR}/CRISPR_26Q1.feather")
print("  -> CRISPR_26Q1.feather written")

# --- XP ---
print("Reading XP_26Q1.rds...")
xp_raw = pyreadr.read_r(f"{DATA_DIR}/XP_26Q1.rds")[None]
# Filter to default model entry (one row per cell line)
xp = xp_raw[xp_raw["IsDefaultEntryForModel"] == "Yes"].copy()
print(f"  XP after default filter: {xp.shape[0]} rows")
gene_cols = [c for c in xp.columns if c not in META_COLS]
xp = xp[["ModelID"] + gene_cols].rename(columns={"ModelID": "row_id"})
xp = xp.reset_index(drop=True)
print(f"  XP: {xp.shape[0]} rows x {xp.shape[1]-1} gene cols")
print(f"  row_id sample: {xp['row_id'].head(3).tolist()}")
xp.to_feather(f"{DATA_DIR}/XP_26Q1.feather")
print("  -> XP_26Q1.feather written")

print("\nDone.")

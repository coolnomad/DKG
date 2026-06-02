"""CLR preprocessing for compositional microbiome OTU tables.

Reads a tab-separated OTU table (taxa x samples), applies centered log-ratio
transformation per sample, then writes two feather files DKG can ingest:

  X.feather  — samples x taxa (all taxa, CLR-transformed), row_id = sample ID
  Y.feather  — samples x 1   (target taxon CLR values),   row_id = sample ID

Usage:
    python clr_preprocess.py \\
        --otu   C:/GitHub/microbiome_studies/14261087/16S_profiles.zip:16S/Taxonomy_16S.txt \\
        --meta  C:/GitHub/microbiome_studies/14261087/Metadata_All.zip:Metadata/metadata.16S.txt \\
        --target g__Akkermansia \\
        --out   C:/GitHub/microbiome_studies/14261087/clr

Zip paths use the colon syntax above; plain .txt/.tsv paths also accepted.

Metadata join:
    The metadata file is assumed to have the same column order as the OTU table
    (positional alignment by column index). The PD/HC label column is written
    into X as a non-numeric annotation column named "diagnosis" — useful for
    phase 11 covariate adjustment but not fed into correlation phases directly.
    Pass --meta-label-col to name it explicitly if the column header differs.

Pseudocount:
    Zeros cannot be log-transformed. A per-sample multiplicative replacement
    is applied before CLR: zeros are replaced with (0.5 * min_nonzero) so
    the relative proportions of non-zero taxa are unchanged.
"""

import argparse
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _open_possibly_zipped(path_str: str) -> io.TextIOWrapper:
    """Open a file that may be inside a zip archive.

    Accepts either a plain path or "archive.zip:inner/path.txt".
    Splits on ".zip:" to avoid tripping on Windows drive letters (C:).
    """
    marker = ".zip:"
    if marker in path_str:
        idx = path_str.index(marker) + len(".zip")
        zip_path = path_str[:idx]
        inner = path_str[idx + 1:]
        zf = zipfile.ZipFile(zip_path)
        return io.TextIOWrapper(zf.open(inner), encoding="utf-8")
    return open(path_str, encoding="utf-8")


def load_otu_table(path_str: str) -> pd.DataFrame:
    """Load OTU table. Returns DataFrame (taxa x samples)."""
    with _open_possibly_zipped(path_str) as fh:
        df = pd.read_csv(fh, sep="\t", index_col=0)
    print(f"OTU table: {df.shape[0]} taxa x {df.shape[1]} samples")
    return df


def load_metadata(path_str: str) -> pd.DataFrame:
    """Load metadata. Returns DataFrame indexed 0..n-1."""
    with _open_possibly_zipped(path_str) as fh:
        raw = fh.read()
    # detect tab-offset header (more values than headers in first data row)
    lines = raw.splitlines()
    header_cols = lines[0].split("\t")
    first_row_cols = lines[1].split("\t")
    if len(first_row_cols) == len(header_cols) + 1:
        # unnamed leading column — insert placeholder
        header_cols = ["_row_num"] + header_cols
    meta = pd.read_csv(
        io.StringIO(raw),
        sep="\t",
        header=0,
        names=header_cols,
    )
    print(f"Metadata: {meta.shape[0]} rows x {meta.shape[1]} cols")
    return meta


# ---------------------------------------------------------------------------
# CLR transformation
# ---------------------------------------------------------------------------

def multiplicative_replacement(mat: np.ndarray) -> np.ndarray:
    """Replace zeros with 0.5 * per-sample min-nonzero value.

    mat shape: (n_samples, n_taxa). Operates row-wise.
    """
    out = mat.copy().astype(np.float64)
    for i in range(out.shape[0]):
        row = out[i]
        nonzero = row[row > 0]
        if len(nonzero) == 0:
            # all-zero sample — leave as-is; will be dropped downstream
            continue
        pseudocount = 0.5 * nonzero.min()
        row[row == 0] = pseudocount
        out[i] = row
    return out


def clr_transform(mat: np.ndarray) -> np.ndarray:
    """Centered log-ratio transformation.

    mat shape: (n_samples, n_taxa), strictly positive after pseudocount.
    Returns same shape in real space.
    """
    log_mat = np.log(mat)
    gm = log_mat.mean(axis=1, keepdims=True)  # geometric mean in log space
    return log_mat - gm


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--otu",    required=True, help="OTU table path (taxa x samples, TSV or zip:inner)")
    parser.add_argument("--meta",   required=True, help="Metadata path (TSV or zip:inner)")
    parser.add_argument("--target", required=True, help="Target taxon name (row in OTU table, e.g. g__Akkermansia)")
    parser.add_argument("--out",    required=True, help="Output directory for X.feather and Y.feather")
    parser.add_argument("--meta-label-col", default="PD", help="Metadata column holding PD/HC label (default: PD)")
    parser.add_argument("--meta-study-col", default="Study_2", help="Metadata column holding cohort/study ID (default: Study_2)")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- load ---
    otu = load_otu_table(args.otu)           # (taxa x samples)
    meta = load_metadata(args.meta)          # (samples x meta_cols)

    # --- validate target ---
    if args.target not in otu.index:
        close = [t for t in otu.index if args.target.lower() in t.lower()]
        msg = f"Target taxon {args.target!r} not found in OTU table."
        if close:
            msg += f" Did you mean one of: {close[:5]}"
        sys.exit(msg)

    # --- align metadata to OTU columns positionally ---
    n_samples = otu.shape[1]
    if len(meta) != n_samples:
        sys.exit(
            f"Metadata rows ({len(meta)}) != OTU columns ({n_samples}). "
            "Cannot positionally align."
        )
    sample_ids = otu.columns.tolist()
    meta = meta.copy()
    meta.index = sample_ids

    # --- drop all-zero samples ---
    col_sums = otu.sum(axis=0)
    zero_samples = col_sums[col_sums == 0].index.tolist()
    if zero_samples:
        print(f"Dropping {len(zero_samples)} all-zero samples")
        otu = otu.drop(columns=zero_samples)
        meta = meta.drop(index=zero_samples)

    # --- transpose to (samples x taxa) for CLR ---
    mat = otu.T.values.astype(np.float64)   # (n_samples, n_taxa)
    taxa_names = otu.index.tolist()
    sample_ids = otu.columns.tolist()

    print(f"Applying multiplicative replacement for zeros...")
    mat = multiplicative_replacement(mat)

    print(f"Applying CLR transform...")
    mat_clr = clr_transform(mat)            # (n_samples, n_taxa)

    # --- split into X and Y ---
    target_idx = taxa_names.index(args.target)
    y_vals = mat_clr[:, target_idx]

    x_taxa = [t for t in taxa_names if t != args.target]
    x_idx  = [i for i, t in enumerate(taxa_names) if t != args.target]
    x_vals = mat_clr[:, x_idx]             # (n_samples, n_taxa - 1)

    # --- build output DataFrames ---
    X_df = pd.DataFrame(x_vals, index=sample_ids, columns=x_taxa)
    X_df.index.name = "row_id"

    Y_df = pd.DataFrame({args.target: y_vals}, index=sample_ids)
    Y_df.index.name = "row_id"

    # annotation columns go to a separate meta file — keeps X purely numeric
    meta_cols: dict[str, list] = {}
    if args.meta_label_col in meta.columns:
        meta_cols["diagnosis"] = meta[args.meta_label_col].tolist()
    if args.meta_study_col in meta.columns:
        meta_cols["study"] = meta[args.meta_study_col].tolist()
    meta_out = pd.DataFrame(meta_cols, index=sample_ids)
    meta_out.index.name = "row_id"

    # --- write ---
    x_path    = out_dir / "X.feather"
    y_path    = out_dir / "Y.feather"
    meta_path = out_dir / "meta.feather"
    X_df.reset_index().to_feather(str(x_path))
    Y_df.reset_index().to_feather(str(y_path))
    meta_out.reset_index().to_feather(str(meta_path))

    print(f"\nX    -> {x_path}  ({X_df.shape[0]} samples x {X_df.shape[1]} taxa, CLR-transformed)")
    print(f"Y    -> {y_path}  ({Y_df.shape[0]} samples x 1 col: {args.target})")
    print(f"meta -> {meta_path}  (diagnosis, study — for covariate adjustment)")
    print("\nDKG usage:")
    print(f"  dkg run --x {x_path} --y {y_path} --y-col {args.target!r} ...")


if __name__ == "__main__":
    main()

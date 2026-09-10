"""
Build CytoSPACE's --n-cells-per-spot-path (-ncpsp) file for each HD sample.

Why this exists: CytoSPACE's default cell-count-per-spot estimation (and its
-mcn/--mean-cell-numbers global-mean flag) is tuned around 55um Visium spots.
At 16um bins, information content per bin varies far more than across 55um
spots (median genes/bin 109-207, but the p10/p90 spread within a sample is
35-780 — see plan Part D), so a single global mean represents that poorly.
This derives a per-bin estimate instead:

    cells_per_bin_i = max(1, round(total_UMI_bin_i / median_total_UMI_per_cell_in_scRNA_ref))

using files already produced by export_cytospace_inputs_hd.py (ST_counts.txt)
and the symlinked scRNA_counts.txt (shared reference, same for every sample).

Output format (verified against the installed CytoSpace package source,
cytospace/common/common.py:read_file + cytospace/cytospace.py:read_data): a
tab-separated file, first column = SpotID (must match ST_coordinates.txt's
SpotID values exactly), second column = an integer cell count. CytoSpace reads
the file with the first column as the row index and takes the first data
column as the count.

MANDATORY CHECKPOINT (plan Part G, #8): before trusting this for a full run,
validate the file against a truncated ~500-bin CytoSPACE run and confirm it
parses without an IndexError about mismatched spot IDs.

NOT YET RUN. Run with: conda run -n squidpy python3 derive_ncells_per_spot.py
(after export_cytospace_inputs_hd.py has produced ST_counts.txt for both
samples; needs the scRNA_counts.txt symlink to already resolve).
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
CYTO_IN = f"{BASE}/cytospace_inputs"
SAMPLES = sys.argv[1:] or ["P18", "P11"]
SC_COUNTS_PATH = f"{CYTO_IN}/scRNA_counts.txt"  # symlinked from reg_visium project


def median_total_umi_per_cell(sc_counts_path, chunksize=2000):
    """Column sums (per-cell total UMI) computed in row chunks to avoid
    loading the full genes x cells matrix into memory at once."""
    col_sums = None
    for chunk in pd.read_csv(sc_counts_path, sep="\t", index_col=0, chunksize=chunksize):
        s = chunk.sum(axis=0)
        col_sums = s if col_sums is None else col_sums.add(s, fill_value=0)
    return float(np.median(col_sums.values))


print(f"Computing median per-cell total UMI from {SC_COUNTS_PATH} (this reads "
      f"the full ~95k-cell reference in chunks — may take a few minutes)...")
median_umi_per_cell = median_total_umi_per_cell(SC_COUNTS_PATH)
print(f"Median total UMI per scRNA cell: {median_umi_per_cell:.2f}")

for sample in SAMPLES:
    st_counts_path = f"{CYTO_IN}/{sample}/ST_counts.txt"
    out_path = f"{CYTO_IN}/{sample}/ST_ncells_per_spot.txt"

    print(f"\n{sample}: reading {st_counts_path}...")
    # genes x bins — sum each column (bin) to get total UMI per bin, in
    # chunks over the gene (row) axis to bound memory
    bin_totals = None
    for chunk in pd.read_csv(st_counts_path, sep="\t", index_col=0, chunksize=2000):
        s = chunk.sum(axis=0)
        bin_totals = s if bin_totals is None else bin_totals.add(s, fill_value=0)

    n_cells = np.maximum(1, np.round(bin_totals / median_umi_per_cell)).astype(int)
    out_df = pd.DataFrame({"n_cells": n_cells})
    out_df.index.name = "SpotID"
    out_df.to_csv(out_path, sep="\t")
    print(f"  {len(out_df)} bins | cells/bin distribution: "
          f"min={n_cells.min()} median={int(np.median(n_cells))} max={n_cells.max()}")
    print(f"  written -> {out_path}")

print("\nDone. Validate on a truncated run before the full CytoSPACE launch (plan Part G, #8).")

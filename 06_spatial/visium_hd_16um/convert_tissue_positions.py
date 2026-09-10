"""
One-off conversion: tissue_positions.parquet -> tissue_positions.csv

Giotto's native Visium HD loaders (createGiottoVisiumHDObject / importVisiumHD)
require the `arrow` R package to read the parquet file, which is not installed
here. Rather than install it, we bypass Giotto's HD-specific codepath entirely:
convert the parquet to a csv once, then feed it into the SAME
createGiottoVisiumObject() call the regular-Visium scripts already use. The
column schema is identical to regular Visium's tissue_positions.csv
(barcode, in_tissue, array_row, array_col, pxl_row_in_fullres,
pxl_col_in_fullres), so no reshaping is needed - this maximises parity between
the 55um and 16um pipelines (see plan Part B).

Run in the `squidpy` conda env (has pandas + pyarrow):
    conda run -n squidpy python3 convert_tissue_positions.py
"""
import os
import sys

import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = sys.argv[1:] or ["P18", "P11"]

for sample in SAMPLES:
    parquet_path = f"{BASE}/sample_{sample}/spatial/tissue_positions.parquet"
    csv_path = f"{BASE}/sample_{sample}/spatial/tissue_positions.csv"

    print(f"{sample}: reading {parquet_path}")
    df = pd.read_parquet(parquet_path)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    df.to_csv(csv_path, index=False)
    print(f"  written -> {csv_path}")

print("\nDone. Verify row counts match the parquet files before proceeding "
      "(see plan Part G, checkpoint 1).")

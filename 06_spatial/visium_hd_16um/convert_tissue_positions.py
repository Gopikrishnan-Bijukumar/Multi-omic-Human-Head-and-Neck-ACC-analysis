"""
Convert tissue_positions.parquet -> tissue_positions.csv for each Visium HD sample.

The csv has the same columns as regular Visium's tissue_positions.csv
(barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres),
so the HD data can be loaded with createGiottoVisiumObject() without the `arrow`
R package. Requires pandas + pyarrow.

    python3 convert_tissue_positions.py [SAMPLE ...]
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

print("\nDone.")

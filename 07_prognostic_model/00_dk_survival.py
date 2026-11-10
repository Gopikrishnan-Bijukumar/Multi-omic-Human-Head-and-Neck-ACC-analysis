"""Build the DK overall-survival table (time in months, event) used by the figure scripts."""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import pandas as pd

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
os.makedirs(OUT, exist_ok=True)

dk_z = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_panel30_z.csv", index_col=0)
surv = pd.read_csv(f"{ROOT}/external_dk/DK_sample_data_table2.csv").set_index("UID")

# patients with expression and usable survival
s = surv.loc[[i for i in dk_z.index if i in surv.index]]
ok = s["SURV_MO"].notna() & s["SURV_CENS"].notna() & (s["SURV_MO"].astype(str) != "NA")
s = s[ok]
time = pd.to_numeric(s["SURV_MO"])
event = pd.to_numeric(s["SURV_CENS"])
print("n with survival:", len(s), " events:", int(event.sum()), " median follow-up:", time.median())

pd.DataFrame({"time": time, "event": event}).to_csv(f"{OUT}/dk_survival.csv")

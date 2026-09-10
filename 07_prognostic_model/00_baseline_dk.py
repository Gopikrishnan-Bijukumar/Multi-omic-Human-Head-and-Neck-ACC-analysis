"""Baseline: apply the existing elastic-net final model to DK, look at survival split.

READ-ONLY on all project inputs. Writes only under work/outputs.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json, numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"

fm = json.load(open(f"{ROOT}/outputs/04_model/final_model.json"))
dk_z = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_panel30_z.csv", index_col=0)
dk_zi = pd.read_csv(f"{ROOT}/outputs/03_harmonized/dk_panel30_z_internal.csv", index_col=0)
surv = pd.read_csv(f"{ROOT}/external_dk/DK_sample_data_table2.csv")
surv = surv.set_index("UID")

print("dk_z", dk_z.shape, "dk_zi", dk_zi.shape)
print(surv[["STATUS", "SURV_CENS", "SURV_MO", "STAGE", "FORM"]].head())

# usable survival
s = surv.loc[[i for i in dk_z.index if i in surv.index]]
ok = s["SURV_MO"].notna() & s["SURV_CENS"].notna() & (s["SURV_MO"].astype(str) != "NA")
print("n with survival:", ok.sum(), "of", len(s))
s = s[ok]
time = pd.to_numeric(s["SURV_MO"])
event = pd.to_numeric(s["SURV_CENS"])
print("events:", int(event.sum()), "median fu:", time.median())

coef = pd.Series(fm["coefficients"])
coef = coef[coef != 0]
print("\nnonzero coefs:\n", coef)

for name, X in [("train_scaler", dk_z), ("dk_internal_scaler", dk_zi)]:
    xs = X.loc[s.index, coef.index]
    score = xs.values @ coef.values + fm["intercept"]
    df = pd.DataFrame({"time": time.values, "event": event.values, "score": score}, index=s.index)
    cph = CoxPHFitter().fit(df[["time", "event", "score"]], "time", "event")
    hi = df["score"] > df["score"].median()
    lr = logrank_test(df.time[hi], df.time[~hi], df.event[hi], df.event[~hi])
    c = cph.concordance_index_
    print(f"\n=== {name} ===  Cox p={cph.summary['p'].iloc[0]:.4g}  HR={np.exp(cph.params_.iloc[0]):.3f}  C={c:.3f}  median-split logrank p={lr.p_value:.4g}")

pd.DataFrame({"time": time, "event": event}).to_csv(f"{OUT}/dk_survival.csv")

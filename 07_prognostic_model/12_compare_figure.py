"""Side-by-side: old elastic-net model vs the new locked signature, on the same DK patients."""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

D = d3.load("all")
T, E = D["T"], D["E"]
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
idx = D["Zd"].index

old = json.load(open(f"{R}/outputs/04_model/final_model.json"))
dkz = pd.read_csv(f"{R}/outputs/03_harmonized/dk_panel30_z.csv", index_col=0).loc[idx]
co = pd.Series(old["coefficients"]); co = co[co != 0]
old_score = dkz[co.index].values @ co.values

fm = json.load(open(f"{OUT}/FINAL_model.json"))
new_score = np.load(f"{OUT}/dk_final_score.npy")

rows = []
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for j, (name, sc) in enumerate([("Previous model\n(elastic net, 7 genes)", old_score),
                                ("New model\n(evidence-informed Bayesian shrinkage)", new_score)]):
    for i, (tag, T_, E_) in enumerate([("5-year OS (primary)", T5, E5), ("full follow-up", T, E)]):
        ax = axes[i, j]
        hi = sc > np.median(sc)
        lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
        cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
        hr = np.exp(cf.params_.iloc[0]); ci = np.exp(cf.confidence_intervals_.iloc[0].values)
        c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
        for m_, lab, col in [(hi, "high risk", "#c0392b"), (~hi, "low risk", "#2471a3")]:
            KaplanMeierFitter().fit(T_[m_], E_[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, color=col)
        ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("overall survival")
        ax.set_title(f"{name}\n{tag}: C={c:.3f}  HR={hr:.2f} [{ci[0]:.2f}-{ci[1]:.2f}]  p={lr.p_value:.3f}", fontsize=9)
        ax.legend(fontsize=7)
        rows.append(dict(model=name.replace("\n", " "), endpoint=tag, C=c, hr=hr,
                         hr_lo=ci[0], hr_hi=ci[1], logrank_p=lr.p_value))
fig.tight_layout(); fig.savefig(f"{OUT}/FINAL_comparison.png", dpi=150)
res = pd.DataFrame(rows)
print(res.to_string(index=False))
res.to_csv(f"{OUT}/FINAL_comparison.csv", index=False)
print("\ncorrelation old vs new score:", round(np.corrcoef(old_score, new_score)[0, 1], 3))

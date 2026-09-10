"""(a) Is DK expression dominated by technical variation? Remove it WITHOUT using survival.
(b) Do independent pre-specified biological axes combine into a stronger score?
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD

D = d3.load("all")
genes, T, E = D["genes"], D["T"], D["E"]
Zd = D["Zd"].values.copy()
clin = D["clin"]
M, mnames, msrc = MOD.build(genes, R)
print("modules", M.shape)


def ev(sc, name):
    sc = np.asarray(sc, float)
    if np.std(sc) < 1e-12: return None
    c = concordance_index_censored(E.astype(bool), T, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "score": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    hrs = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "g": hi.astype(int)}), "time", "event").params_.iloc[0])
    return dict(score=name, C=c, hr_per_sd=np.exp(cph.params_.iloc[0]), cox_p=cph.summary["p"].iloc[0],
                logrank_p=lr.p_value, hr_split=hrs)


# ---------- technical covariates ----------
tech = clin[["TOTAL", "COUNTS", "OFF_TARG", "AMBIG", "LOWQ", "UNALIGN"]].apply(pd.to_numeric, errors="coerce")
tech["frac_counts"] = tech["COUNTS"] / tech["TOTAL"]
tech["frac_offt"] = tech["OFF_TARG"] / tech["TOTAL"]
tech["frac_lowq"] = tech["LOWQ"] / tech["TOTAL"]
techm = tech[["frac_counts", "frac_offt", "frac_lowq", "TOTAL"]].copy()
techm["logtotal"] = np.log10(techm.pop("TOTAL"))
techz = ((techm - techm.mean()) / techm.std()).fillna(0).values

U, S, Vt = np.linalg.svd(Zd - Zd.mean(0), full_matrices=False)
pcs = U[:, :10] * S[:10]
print("\nvariance explained:", (S[:10] ** 2 / (S ** 2).sum()).round(3))
print("PC vs technical covariate |r| (rows PC1-10, cols frac_counts/frac_offt/frac_lowq/logtotal):")
cor = np.array([[np.corrcoef(pcs[:, i], techz[:, j])[0, 1] for j in range(techz.shape[1])] for i in range(10)])
print(np.round(cor, 2))
print("PC vs survival (cox p):")
for i in range(6):
    r = ev(pcs[:, i], f"PC{i+1}")
    print(f"  PC{i+1}: C={r['C']:.3f} p={r['cox_p']:.3g}  (tech |r|max={np.abs(cor[i]).max():.2f})")

# genes' loading on technical PCs -> remove PCs that are technical (|r|>0.5 with any covariate)
tech_pcs = [i for i in range(10) if np.abs(cor[i]).max() > 0.5]
print("\ntechnical PCs removed:", tech_pcs)
Zd_c = Zd - Zd.mean(0)
if tech_pcs:
    for i in tech_pcs:
        v = Vt[i]
        Zd_c = Zd_c - np.outer(Zd_c @ v, v)
Zd_dn = stats.zscore(Zd_c, axis=0)

# ---------- axis scores, raw vs denoised ----------
axes = {"hypoxia": "hypoxia", "glycolysis": "glycolysis", "MYC_translation": "MYC_translation",
        "cell_cycle_G2M": "cell_cycle_G2M", "cell_cycle_S": "cell_cycle_S", "E2F_targets": "E2F_targets",
        "EMT": "EMT", "evidence_poor": "evidence_poor", "evidence_good": "evidence_good",
        "wgcna_blue": "wgcna_blue", "wgcna_turquoise": "wgcna_turquoise",
        "sc_Actively_Dividing_cells_up_poor": "sc_Actively_Dividing_cells_up_poor"}
rows = []
sc_raw, sc_dn = {}, {}
for a in axes:
    if a not in mnames: print("missing module", a); continue
    k = mnames.index(a)
    sc_raw[a] = Zd @ M[:, k]; sc_dn[a] = Zd_dn @ M[:, k]
    for tag, sdict in [("raw", sc_raw), ("denoised", sc_dn)]:
        r = ev(sdict[a], f"{a}[{tag}]")
        if r: rows.append(r)

print("\n--- correlation between key axes (raw) ---")
key = [a for a in ["hypoxia", "MYC_translation", "evidence_poor", "cell_cycle_G2M", "wgcna_blue"] if a in sc_raw]
print(pd.DataFrame({a: sc_raw[a] for a in key}).corr().round(2).to_string())

# ---------- fixed equal-weight combinations (no fitting) ----------
def z(x): return stats.zscore(np.asarray(x, float))
combos = {
 "hyp+evid":        lambda s: z(s["hypoxia"]) + z(s["evidence_poor"]) - z(s["evidence_good"]),
 "hyp+myc":         lambda s: z(s["hypoxia"]) + z(s["MYC_translation"]),
 "hyp+myc+evid":    lambda s: z(s["hypoxia"]) + z(s["MYC_translation"]) + z(s["evidence_poor"]) - z(s["evidence_good"]),
 "hyp+prolif":      lambda s: z(s["hypoxia"]) + z(s["cell_cycle_G2M"]),
 "hyp+myc+prolif":  lambda s: z(s["hypoxia"]) + z(s["MYC_translation"]) + z(s["cell_cycle_G2M"]),
 "evid_only":       lambda s: z(s["evidence_poor"]) - z(s["evidence_good"]),
}
for name, fn in combos.items():
    for tag, sdict in [("raw", sc_raw), ("denoised", sc_dn)]:
        try: r = ev(fn(sdict), f"COMBO {name}[{tag}]")
        except KeyError: continue
        if r: rows.append(r)

res = pd.DataFrame(rows).sort_values("C", ascending=False)
pd.set_option("display.width", 200)
print("\n", res.to_string(index=False))
res.to_csv(f"{ACC_DATA_ROOT}/model_outputs/denoise_combine.csv", index=False)
np.save(f"{ACC_DATA_ROOT}/model_outputs/Zd_denoised.npy", Zd_dn)

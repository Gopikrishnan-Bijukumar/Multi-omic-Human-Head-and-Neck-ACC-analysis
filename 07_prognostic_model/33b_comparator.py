"""Stage 33b - penalized-regression comparator, done properly (review item #1).

Stage 33 part E used a single fixed C = 0.1 for every penalty. That is degenerate for L1
on 12,700 features (all coefficients shrink to zero, LOO AUC 0.500) and arbitrary for the
others, so it cannot answer the review's question "why variational inference rather than a
simpler penalized model?".

This stage answers it under the same protocol the Bayesian model was held to:

  1. choose the penalty strength by leave-one-out AUC INSIDE the 20 discovery patients;
  2. refit on all 20;
  3. evaluate once on DK, both as the dense linear score and as a top-4 sign panel.

DK survival is used only in step 3, exactly as for the locked panel.

Supersedes the `comparator` block of outputs_vi/VI_DIAGNOSTICS.json, which is rewritten
in place with these results (that file is stage-33 output, not a locked artefact).
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util, os, time
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sksurv.metrics import concordance_index_censored
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_vi"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
keep = (dk_raw.mean(axis=0).values > 1.0) & ((dk_raw > 0).mean(axis=0).values > 0.6) & (tr_raw.mean(axis=0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gidx = {g: i for i, g in enumerate(genes)}
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
print(f"genes={len(genes)} JSE n={len(y)} DK n={len(T)} events5y={E5.sum()}", flush=True)


def mk(kind, C):
    if kind == "l1":
        return LogisticRegression(penalty="l1", C=C, solver="liblinear", max_iter=20000)
    if kind == "l2":
        return LogisticRegression(penalty="l2", C=C, solver="lbfgs", max_iter=20000)
    return LogisticRegression(penalty="elasticnet", C=C, l1_ratio=0.5, solver="saga", max_iter=8000)


def dk_c(sc):
    return float(concordance_index_censored(E5.astype(bool), T5, sc)[0])


GRID = {"l1": [0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
        "l2": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
        "enet": [0.1, 0.3, 1.0, 3.0]}
rows, t0 = [], time.time()
for kind, Cs in GRID.items():
    for C in Cs:
        pred = np.zeros(len(y))
        for i in range(len(y)):
            tr = np.array([j for j in range(len(y)) if j != i])
            clf = mk(kind, C).fit(Zt[tr], y[tr])
            pred[i] = clf.decision_function(Zt[i:i + 1])[0]
        loo = float(roc_auc_score(y, pred))
        clf = mk(kind, C).fit(Zt, y)
        co = clf.coef_.ravel()
        nz = int((co != 0).sum())
        dense = dk_c(Zd @ co) if nz else 0.5
        if nz >= 4:
            ii = np.argsort(-np.abs(co))[:4]
            pan = [genes[i] for i in ii]
            top4 = dk_c((Zd[:, ii] * np.sign(co[ii])).mean(1))
        else:
            pan, top4 = [], 0.5
        rows.append({"penalty": kind, "C": C, "loo_auc": loo, "n_nonzero": nz,
                     "dk_c_dense": dense, "dk_c_top4": top4, "top4": ";".join(pan),
                     "n_locked_in_top4": len(set(pan) & set(PANEL4))})
        print(f"  {kind:5s} C={C:<6g} LOO {loo:.3f}  nz {nz:6d}  DK C dense {dense:.3f}  top4 {top4:.3f}", flush=True)
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/comparator_path.csv", index=False)

# selection rule: best discovery LOO AUC, ties broken toward the sparser model
best = df.sort_values(["loo_auc", "n_nonzero"], ascending=[False, True]).iloc[0]
print(f"\nselected by discovery LOO only: {best.penalty} C={best.C} "
      f"(LOO {best.loo_auc:.3f}, nz {best.n_nonzero})", flush=True)
print(f"  -> DK 5-yr C  dense {best.dk_c_dense:.3f} | top-4 {best.dk_c_top4:.3f} | "
      f"locked panel 0.671", flush=True)

comp = {
    "protocol": ("penalty strength chosen by leave-one-out AUC inside the 20 discovery patients, "
                 "refit on all 20, evaluated once on DK - the same protocol as the locked panel"),
    "path": df.to_dict("records"),
    "selected": {"penalty": best.penalty, "C": float(best.C), "loo_auc": float(best.loo_auc),
                 "n_nonzero": int(best.n_nonzero), "dk_c_dense": float(best.dk_c_dense),
                 "dk_c_top4": float(best.dk_c_top4), "top4": best.top4,
                 "n_locked_in_top4": int(best.n_locked_in_top4)},
    "best_dk_c_dense_anywhere_on_path": float(df.dk_c_dense.max()),
    "best_dk_c_top4_anywhere_on_path": float(df.dk_c_top4.max()),
    "locked_panel_dk_c": 0.67149,
    "bayesian_loo_auc_reported": 0.96,
    "why_not_mcmc": ("The model carries 12,700 gene-level and 52 module-level latents plus "
                     "hierarchical scales - >25,000 parameters fitted to 20 binary outcomes. "
                     "MCMC was not attempted. The 3000-step run in stage 33 part C is the "
                     "convergence reference. This is a stated limitation, not a claim that the "
                     "variational posterior is accurate."),
    "supersedes": "the `comparator` block of stage 33 part E, which used a single fixed C=0.1"}

res = json.load(open(f"{OUT}/VI_DIAGNOSTICS.json"))
res["comparator_superseded_by_33b"] = res.pop("comparator")
res["comparator"] = comp
json.dump(res, open(f"{OUT}/VI_DIAGNOSTICS.json", "w"), indent=2)
json.dump(comp, open(f"{OUT}/COMPARATOR.json", "w"), indent=2)
print(f"\nwrote {OUT}/COMPARATOR.json and updated VI_DIAGNOSTICS.json ({time.time()-t0:.0f}s)", flush=True)

"""Repeated cross-validated evaluation of the multi-level model on DK.

Every model is refitted inside every fold; DK survival never leaves the training folds.
The module definitions, the evidence weights and the JSE cohort are DK-independent, so
they enter each fold in full.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, time, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, prep_cox, to_t

N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 5
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 1500

D = d3.load("all")
genes_all = D["genes"]
Zd_all, Zt_all = D["Zd"], D["Zt"]
T, E = D["T"], D["E"]

M_all, mod_names, mod_src = MOD.build(genes_all, R)
keep = np.abs(M_all).sum(1) > 0                     # genes that belong to >=1 module
# plus the evidence candidate genes so gene-level deviations can act on them
ev = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
keep |= np.isin(genes_all, ev.index.values)
genes = [g for g, k in zip(genes_all, keep) if k]
idx = np.where(keep)[0]
M = M_all[idx]
Zd, Zt = Zd_all.values[:, idx], Zt_all.values[:, idx]
y = D["y"]
s_ev, d_ev = d3.evidence_prior(ev.reindex(genes))
srcs = sorted(set(mod_src)); src_ids = np.array([srcs.index(s) for s in mod_src])
P, K = M.shape
print(f"genes={P} modules={K} sources={srcs} DK n={len(T)} events={E.sum()}")
mod_scores_dk = Zd @ M                              # (n, K) module scores

tstat = stats.ttest_ind(Zt[y == 1], Zt[y == 0], axis=0).statistic
tp = stats.ttest_ind(Zt[y == 1], Zt[y == 0], axis=0).pvalue
Xj, yj = to_t(Zt), to_t(y)


def cidx(sc): return concordance_index_censored(E.astype(bool), T, sc)[0]


def fit_ml(tr, **kw):
    Xd, evv = prep_cox(Zd[tr], T[tr], E[tr])
    m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, seed=kw.pop("seed", 0), **kw)
    m.fit(Xj=Xj, yj=yj, Xd=Xd, event=evv, steps=STEPS)
    return m.posterior(300)


def mdl_full(tr, te):    return Zd[te] @ fit_ml(tr)["beta"]
def mdl_modonly(tr, te): return Zd[te] @ fit_ml(tr, use_gene_level=False)["beta"]
def mdl_geneonly(tr, te):return Zd[te] @ fit_ml(tr, use_module_level=False)["beta"]
def mdl_nojse(tr, te):   return Zd[te] @ fit_ml(tr, use_jse=False)["beta"]
def mdl_flat(tr, te):
    return Zd[te] @ fit_ml(tr, s_evidence_override=None, use_dir_prior=False)["beta"] if False else \
           Zd[te] @ _flat(tr)


def _flat(tr):
    Xd, evv = prep_cox(Zd[tr], T[tr], E[tr])
    m = MultiLevelModel(M, src_ids, len(srcs), np.ones(P, np.float32), np.zeros(P, np.float32),
                        use_dir_prior=False)
    m.fit(Xj=Xj, yj=yj, Xd=Xd, event=evv, steps=STEPS)
    return m.posterior(300)["beta"]


def mdl_hypoxia(tr, te):
    k = mod_names.index("hypoxia")
    return mod_scores_dk[te, k]


def mdl_modridge(tr, te):
    """Classical comparator: ridge Cox on the K module scores only."""
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv
    Xtr = stats.zscore(mod_scores_dk[tr], axis=0)
    mu, sd = mod_scores_dk[tr].mean(0), mod_scores_dk[tr].std(0) + 1e-9
    best, bestc = None, -1
    for a in [3.0, 1.0, 0.3, 0.1, 0.03]:
        try:
            mm = CoxnetSurvivalAnalysis(l1_ratio=1e-8, alphas=[a], max_iter=50000)
            mm.fit(Xtr, Surv.from_arrays(E[tr].astype(bool), T[tr]))
            c = concordance_index_censored(E[tr].astype(bool), T[tr], mm.predict(Xtr))[0]
        except Exception: continue
        if c > bestc: best, bestc = mm, c
    if best is None: return np.zeros(len(te))
    return best.predict((mod_scores_dk[te] - mu) / sd)


def mdl_jse_sign(tr, te, k=100):
    o = np.argsort(tp)[:k]
    b = np.zeros(P); b[o] = np.sign(tstat[o])
    return Zd[te] @ b


def mdl_evdir(tr, te):
    return Zd[te] @ (d_ev * s_ev)


MODELS = {
    "ML_bayes_full":       mdl_full,
    "ML_bayes_moduleonly": mdl_modonly,
    "ML_bayes_geneonly":   mdl_geneonly,
    "ML_bayes_noJSE":      mdl_nojse,
    "ML_bayes_flatprior":  mdl_flat,
    "module_ridge_cox":    mdl_modridge,
    "hypoxia_fixed":       mdl_hypoxia,
    "jse_sign_top100":     mdl_jse_sign,
    "evidence_direction":  mdl_evdir,
}

oof = {k: np.zeros((N_REP, len(T))) for k in MODELS}
percv = {k: [] for k in MODELS}
t0 = time.time()
for rep in range(N_REP):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=100 + rep).split(Zd, E):
        for name, fn in MODELS.items():
            try: sc = np.asarray(fn(tr, te), float)
            except Exception as ex:
                print("FAIL", name, ex); sc = np.zeros(len(te))
            oof[name][rep, te] = stats.zscore(sc) if np.std(sc) > 1e-12 else 0.0
    for k in MODELS: percv[k].append(cidx(oof[k][rep]))
    print(f"rep {rep} {time.time()-t0:.0f}s " + " ".join(f"{k}={percv[k][-1]:.3f}" for k in MODELS))

rows = []
for name in MODELS:
    m = oof[name].mean(0)
    hi = m > np.median(m)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "score": stats.zscore(m)}), "time", "event")
    d2 = pd.DataFrame({"time": T, "event": E, "grp": hi.astype(int)})
    hrs = np.exp(CoxPHFitter().fit(d2, "time", "event").params_.iloc[0])
    rows.append(dict(model=name, c_mean=np.mean(percv[name]), c_sd=np.std(percv[name]), c_pooled=cidx(m),
                     hr_per_sd=np.exp(cph.params_.iloc[0]), cox_p=cph.summary["p"].iloc[0],
                     logrank_p=lr.p_value, hr_median_split=hrs))
res = pd.DataFrame(rows).sort_values("c_pooled", ascending=False)
pd.set_option("display.width", 220)
print(res.to_string(index=False))
res.to_csv(f"{ACC_DATA_ROOT}/model_outputs/ml_cv_comparison.csv", index=False)
np.savez(f"{ACC_DATA_ROOT}/model_outputs/ml_cv_oof.npz", T=T, E=E, **oof)

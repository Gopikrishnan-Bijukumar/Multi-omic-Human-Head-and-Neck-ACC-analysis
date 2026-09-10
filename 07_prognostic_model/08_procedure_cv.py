"""Cross-validate the whole PROCEDURE, not just the model.

The concern with "hypoxia + MYC scored best on DK" is that the choice of modules was
informed by DK survival. The fix: define the module choice as an algorithm, and put that
algorithm inside the cross-validation loop. Modules are chosen from the DK-independent
library using only the training folds' survival; the held-out fold never contributes.
An out-of-fold C-index from this loop is an honest estimate of the whole pipeline,
selection included.
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

N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 10
D = d3.load("all")
genes, T, E, y = D["genes"], D["T"], D["E"], D["y"]
Zd, Zt = D["Zd"].values, D["Zt"].values
M, mnames, msrc = MOD.build(genes, R)
srcs = sorted(set(msrc)); src_ids = np.array([srcs.index(s) for s in msrc])
SC = Zd @ M                     # (n, K) DK module scores  -- unsupervised
SCt = Zt @ M                    # JSE module scores
K = M.shape[1]
print(f"DK n={len(T)} events={E.sum()}  modules K={K}")

ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
s_ev, d_ev = d3.evidence_prior(ev_tab.reindex(genes))
Xj, yj = to_t(Zt), to_t(y)


def cox_z(X, T_, E_):
    o = np.argsort(T_); Xs = X[o]; Es = E_[o]
    Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
    n = Xs.shape[0]
    cs = np.cumsum(Xs[::-1], 0)[::-1]; cs2 = np.cumsum((Xs ** 2)[::-1], 0)[::-1]
    nr = np.arange(n, 0, -1)[:, None]
    mr = cs / nr; vr = cs2 / nr - mr ** 2
    m = Es == 1
    return (Xs[m] - mr[m]).sum(0) / np.sqrt(vr[m].sum(0) + 1e-12)


# JSE module-level evidence: t statistic of each module score, Poor vs Good (DK-independent)
tmod = stats.ttest_ind(SCt[y == 1], SCt[y == 0], axis=0).statistic


def make_topk(k, pool=None, use_jse_prior=False, alpha=0.0):
    idx_pool = np.arange(K) if pool is None else np.array([i for i in range(K) if msrc[i] in pool])
    def f(tr, te):
        z = cox_z(SC[tr][:, idx_pool], T[tr], E[tr])
        crit = np.abs(z) + alpha * np.abs(tmod[idx_pool]) if use_jse_prior else np.abs(z)
        top = idx_pool[np.argsort(-crit)[:k]]
        zz = cox_z(SC[tr][:, top], T[tr], E[tr])
        mu, sd = SC[tr][:, top].mean(0), SC[tr][:, top].std(0) + 1e-9
        return (((SC[te][:, top] - mu) / sd) * np.sign(zz)).mean(1)
    return f


def make_shrunk(pool=None, kappa=2.0):
    """Soft-thresholded module weights: w_k = sign(z)*max(|z|-kappa,0)."""
    idx_pool = np.arange(K) if pool is None else np.array([i for i in range(K) if msrc[i] in pool])
    def f(tr, te):
        z = cox_z(SC[tr][:, idx_pool], T[tr], E[tr])
        w = np.sign(z) * np.maximum(np.abs(z) - kappa, 0)
        if np.all(w == 0): w = np.sign(z) * (np.abs(z) == np.abs(z).max())
        mu, sd = SC[tr][:, idx_pool].mean(0), SC[tr][:, idx_pool].std(0) + 1e-9
        return ((SC[te][:, idx_pool] - mu) / sd) @ w / (np.abs(w).sum() + 1e-9)
    return f


def make_bayes_mod(omega0, use_jse=True, tau0=1e-6):
    def f(tr, te):
        Xd, evv = prep_cox(Zd[tr], T[tr], E[tr])
        m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, tau0=tau0, omega0=omega0,
                            use_gene_level=False, use_jse=use_jse)
        m.fit(Xj=Xj, yj=yj, Xd=Xd, event=evv, steps=1200)
        return Zd[te] @ m.posterior(300)["beta"]
    return f


def fixed(name_list):
    ii = [mnames.index(n) for n in name_list]
    def f(tr, te):
        mu, sd = SC[tr][:, ii].mean(0), SC[tr][:, ii].std(0) + 1e-9
        return ((SC[te][:, ii] - mu) / sd).mean(1)
    return f


PROC = {
    "topk1_all":        make_topk(1),
    "topk2_all":        make_topk(2),
    "topk3_all":        make_topk(3),
    "topk5_all":        make_topk(5),
    "topk2_curated":    make_topk(2, pool={"curated"}),
    "topk3_curated":    make_topk(3, pool={"curated"}),
    "topk2_jseprior":   make_topk(2, use_jse_prior=True, alpha=0.5),
    "topk3_jseprior":   make_topk(3, use_jse_prior=True, alpha=0.5),
    "shrunk_k2.0":      make_shrunk(kappa=2.0),
    "shrunk_k1.5":      make_shrunk(kappa=1.5),
    "bayes_mod_w0.5":   make_bayes_mod(0.5),
    "bayes_mod_w0.2":   make_bayes_mod(0.2),
    "FIXED_hyp_myc":    fixed(["hypoxia", "MYC_translation"]),
    "FIXED_hypoxia":    fixed(["hypoxia"]),
    "FIXED_evidence":   fixed(["evidence_poor"]),
}

oof = {k: np.zeros((N_REP, len(T))) for k in PROC}
per = {k: [] for k in PROC}
picked = {k: [] for k in PROC}
t0 = time.time()
for rep in range(N_REP):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=200 + rep).split(Zd, E):
        for name, fn in PROC.items():
            try: sc = np.asarray(fn(tr, te), float)
            except Exception as ex:
                print("FAIL", name, ex); sc = np.zeros(len(te))
            oof[name][rep, te] = stats.zscore(sc) if np.std(sc) > 1e-12 else 0.0
    for k in PROC:
        per[k].append(concordance_index_censored(E.astype(bool), T, oof[k][rep])[0])
    print(f"rep {rep} {time.time()-t0:.0f}s " + " ".join(f"{k}={per[k][-1]:.3f}" for k in list(PROC)[:8]))

rows = []
for name in PROC:
    m = oof[name].mean(0)
    hi = m > np.median(m)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "score": stats.zscore(m)}), "time", "event")
    hrs = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "g": hi.astype(int)}), "time", "event").params_.iloc[0])
    rows.append(dict(procedure=name, c_mean=np.mean(per[name]), c_sd=np.std(per[name]),
                     c_pooled=concordance_index_censored(E.astype(bool), T, m)[0],
                     hr_per_sd=np.exp(cph.params_.iloc[0]), cox_p=cph.summary["p"].iloc[0],
                     logrank_p=lr.p_value, hr_median_split=hrs))
res = pd.DataFrame(rows).sort_values("c_mean", ascending=False)
pd.set_option("display.width", 220)
print(res.to_string(index=False))
res.to_csv(f"{ACC_DATA_ROOT}/model_outputs/procedure_cv.csv", index=False)
np.savez(f"{ACC_DATA_ROOT}/model_outputs/procedure_oof.npz", T=T, E=E, **oof)

# which modules does the selection actually pick, across all folds?
print("\nmodules selected by topk2_all across folds:")
cnt = {}
for rep in range(N_REP):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=200 + rep).split(Zd, E):
        z = cox_z(SC[tr], T[tr], E[tr])
        for i in np.argsort(-np.abs(z))[:2]:
            cnt[mnames[i]] = cnt.get(mnames[i], 0) + 1
sel = pd.Series(cnt).sort_values(ascending=False) / (N_REP * 5)
print(sel.head(12).to_string())
sel.to_csv(f"{ACC_DATA_ROOT}/model_outputs/module_selection_freq.csv")

"""Honest repeated cross-validated comparison of prognostic models on the DK cohort.

Every model is re-fitted inside every fold. The JSE cohort (20 patients) and all
discovery evidence (WGCNA / single-cell DE / bulk DE) are DK-independent, so they
are used in full inside each fold; only DK survival is folded.

Writes only work/outputs.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, time, json, importlib.util
import numpy as np, pandas as pd
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
from bayes import JointShrinkageModel, prep_cox, to_t
from sklearn.model_selection import StratifiedKFold
from sksurv.metrics import concordance_index_censored
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import warnings; warnings.filterwarnings("ignore")

OUT = f"{ACC_DATA_ROOT}/model_outputs"
N_REPEATS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
POOL = sys.argv[2] if len(sys.argv) > 2 else "evidence"

D = d3.load(POOL)
genes, Zd, Zt, y, T, E = D["genes"], D["Zd"].values, D["Zt"].values, D["y"], D["T"], D["E"]
s_ev, d_ev = d3.evidence_prior(D["ev"])
P = len(genes)
print(f"pool={POOL} genes={P} DK n={len(T)} events={E.sum()} JSE n={len(y)}")

# JSE per-gene t statistic (discovery only)
tstat = stats.ttest_ind(Zt[y == 1], Zt[y == 0], axis=0).statistic
tp = stats.ttest_ind(Zt[y == 1], Zt[y == 0], axis=0).pvalue


def cidx(sc, T_, E_):
    return concordance_index_censored(E_.astype(bool), T_, sc)[0]


# ---------------- model zoo: each returns a score for Xte ----------------
def m_bayes(tr_idx, te_idx, **kw):
    Xd, ev_, _ = prep_cox(Zd[tr_idx], T[tr_idx], E[tr_idx])
    m = JointShrinkageModel(P, kw.get("s", s_ev), kw.get("d", d_ev), tau0=kw.get("tau0", 0.05),
                            use_jse=kw.get("use_jse", True), use_dk=kw.get("use_dk", True),
                            use_dir_prior=kw.get("use_dir", True), seed=kw.get("seed", 0))
    m.fit(Xj=to_t(Zt), yj=to_t(y), Xd=Xd, event=ev_, steps=kw.get("steps", 1200))
    b, _, _ = m.beta_mean()
    return Zd[te_idx] @ b, b


def m_jse_sign(tr_idx, te_idx, k=100):
    o = np.argsort(tp)[:k]
    b = np.zeros(P); b[o] = np.sign(tstat[o])
    return Zd[te_idx] @ b, b


def m_evidence_dir(tr_idx, te_idx):
    b = d_ev * s_ev
    return Zd[te_idx] @ b, b


def m_coxnet(tr_idx, te_idx, l1=1.0):
    ytr = Surv.from_arrays(E[tr_idx].astype(bool), T[tr_idx])
    mdl = CoxnetSurvivalAnalysis(l1_ratio=l1 if l1 > 0 else 1e-8, alpha_min_ratio=0.05, n_alphas=30,
                                 max_iter=20000, normalize=False)
    Xtr = Zd[tr_idx]
    mdl.fit(Xtr, ytr)
    # inner 4-fold CV over the alpha path
    alphas = mdl.alphas_
    inner = StratifiedKFold(4, shuffle=True, random_state=1).split(Xtr, E[tr_idx])
    sc = np.zeros(len(alphas)); cnt = np.zeros(len(alphas))
    for itr, ite in inner:
        if E[tr_idx][ite].sum() < 2: continue
        mi = CoxnetSurvivalAnalysis(l1_ratio=l1 if l1 > 0 else 1e-8, alphas=alphas, max_iter=20000, normalize=False)
        try: mi.fit(Xtr[itr], Surv.from_arrays(E[tr_idx][itr].astype(bool), T[tr_idx][itr]))
        except Exception: continue
        for j, a in enumerate(alphas):
            try:
                p = mi.predict(Xtr[ite], alpha=a)
                if np.std(p) < 1e-12: continue
                sc[j] += cidx(p, T[tr_idx][ite], E[tr_idx][ite]); cnt[j] += 1
            except Exception: pass
    best = alphas[np.argmax(np.where(cnt > 0, sc / np.maximum(cnt, 1), 0))]
    b = mdl.coef_[:, list(alphas).index(best)]
    return Zd[te_idx] @ b, b


def m_rsf(tr_idx, te_idx, k=200):
    o = np.argsort(tp)[:k]                       # DK-independent pre-filter
    mdl = RandomSurvivalForest(n_estimators=300, min_samples_leaf=3, max_features="sqrt",
                               random_state=0, n_jobs=8)
    mdl.fit(Zd[tr_idx][:, o], Surv.from_arrays(E[tr_idx].astype(bool), T[tr_idx]))
    b = np.zeros(P)
    return mdl.predict(Zd[te_idx][:, o]), b


def m_univ_dk(tr_idx, te_idx, k=50):
    """Pure DK-driven univariate screen, refit in fold (comparator without evidence)."""
    zz = _cox_z(Zd[tr_idx], T[tr_idx], E[tr_idx])
    o = np.argsort(-np.abs(zz))[:k]
    b = np.zeros(P); b[o] = np.sign(zz[o])
    return Zd[te_idx] @ b, b


def m_hybrid(tr_idx, te_idx, k=50):
    """Evidence-informed screen: rank by JSE evidence, sign from combined evidence+fold DK."""
    zz = _cox_z(Zd[tr_idx], T[tr_idx], E[tr_idx])
    comb = np.abs(tstat) * s_ev
    o = np.argsort(-comb)[:k]
    b = np.zeros(P); b[o] = np.sign(tstat[o]) * 0.5 + np.sign(zz[o]) * 0.5
    return Zd[te_idx] @ b, b


def _cox_z(X, T_, E_):
    o = np.argsort(T_); Xs = X[o]; Es = E_[o]
    n = Xs.shape[0]
    Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
    csum = np.cumsum(Xs[::-1], 0)[::-1]; csum2 = np.cumsum((Xs ** 2)[::-1], 0)[::-1]
    nr = np.arange(n, 0, -1)[:, None]
    mr = csum / nr; vr = csum2 / nr - mr ** 2
    ev = Es == 1
    U = (Xs[ev] - mr[ev]).sum(0); V = vr[ev].sum(0)
    return U / np.sqrt(V + 1e-12)


MODELS = {
    "bayes_joint_evid":      lambda a, b: m_bayes(a, b),
    "bayes_dkonly_evid":     lambda a, b: m_bayes(a, b, use_jse=False),
    "bayes_joint_flat":      lambda a, b: m_bayes(a, b, s=np.ones(P, np.float32), d=np.zeros(P, np.float32), use_dir=False),
    "bayes_jseonly_evid":    lambda a, b: m_bayes(a, b, use_dk=False),
    "jse_sign_top100":       lambda a, b: m_jse_sign(a, b, 100),
    "evidence_direction":    m_evidence_dir,
    "coxnet_lasso":          lambda a, b: m_coxnet(a, b, 1.0),
    "coxnet_ridge":          lambda a, b: m_coxnet(a, b, 0.01),
    "rsf_top200":            lambda a, b: m_rsf(a, b, 200),
    "dk_univ_top50":         lambda a, b: m_univ_dk(a, b, 50),
    "hybrid_evid_dk_top50":  lambda a, b: m_hybrid(a, b, 50),
}

results = {k: [] for k in MODELS}
oof = {k: np.zeros((N_REPEATS, len(T))) for k in MODELS}
t0 = time.time()
for rep in range(N_REPEATS):
    skf = StratifiedKFold(5, shuffle=True, random_state=100 + rep)
    for tr_idx, te_idx in skf.split(Zd, E):
        for name, fn in MODELS.items():
            try:
                sc, _ = fn(tr_idx, te_idx)
            except Exception as ex:
                sc = np.zeros(len(te_idx)); print("FAIL", name, ex)
            oof[name][rep, te_idx] = stats.zscore(sc) if np.std(sc) > 1e-12 else 0.0
    for name in MODELS:
        results[name].append(cidx(oof[name][rep], T, E))
    print(f"repeat {rep} done  {time.time()-t0:.0f}s  " +
          "  ".join(f"{k}={results[k][-1]:.3f}" for k in MODELS))

rows = []
for name in MODELS:
    m = oof[name].mean(0)
    c = cidx(m, T, E)
    hi = m > np.median(m)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    df = pd.DataFrame({"time": T, "event": E, "score": stats.zscore(m)})
    cph = CoxPHFitter().fit(df, "time", "event")
    rows.append(dict(model=name, c_mean_rep=np.mean(results[name]), c_sd=np.std(results[name]),
                     c_pooled=c, hr_per_sd=np.exp(cph.params_.iloc[0]), cox_p=cph.summary["p"].iloc[0],
                     logrank_p=lr.p_value))
res = pd.DataFrame(rows).sort_values("c_pooled", ascending=False)
print(res.to_string(index=False))
res.to_csv(f"{OUT}/cv_model_comparison_{POOL}.csv", index=False)
np.savez(f"{OUT}/cv_oof_{POOL}.npz", **{k: v for k, v in oof.items()}, T=T, E=E)

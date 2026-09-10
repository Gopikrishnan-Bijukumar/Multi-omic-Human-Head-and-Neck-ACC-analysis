"""FINAL: tune priors on the discovery cohort only, lock the model, validate once on DK.

Protocol
--------
1. The evidence-informed Bayesian multi-level shrinkage model is fitted to the JSE
   discovery cohort ALONE (20 patients, Poor vs Good). Its prior hyper-parameters are
   chosen by leave-one-out cross-validation *within JSE*. DK is not touched.
2. The winning configuration is refitted on all 20 JSE patients -> locked gene weights.
3. Those weights are applied to DK once. Primary endpoint: 5-year overall survival
   (administratively censored at 60 months) - the endpoint the discovery labels encode.
4. Nulls: 2000 random gene sets, and 2000 permutations of DK survival.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util, time
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, to_t

rng = np.random.default_rng(6)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]

# unsupervised expression filter (pre-specified: DK depth + discovery detection)
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
srcs = sorted(set(msrc)); src_ids = np.array([srcs.index(s) for s in msrc])
ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
s_ev, d_ev = d3.evidence_prior(ev_tab.reindex(genes))
print(f"genes={len(genes)} modules={M.shape[1]} JSE n={len(y)} DK n={len(T)} events={E.sum()}")

STEPS = 900


def fit_beta(idx, tau0, omega0, gene_level, module_level, seed=0):
    m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, tau0=tau0, omega0=omega0,
                        use_gene_level=gene_level, use_module_level=module_level,
                        use_dk=False, use_jse=True, seed=seed)
    m.fit(Xj=to_t(Zt[idx]), yj=to_t(y[idx]), steps=STEPS)
    return m.posterior(300)


# ---------------- 1. hyper-parameter selection by LOO inside JSE ----------------
GRID = [dict(tau0=t, omega0=o, gene_level=g, module_level=mo)
        for t in [0.02, 0.1] for o in [0.1, 0.5] for g, mo in [(True, True), (False, True), (True, False)]]
loo_rows = []
t0 = time.time()
for cfg in GRID:
    pred = np.zeros(len(y))
    for i in range(len(y)):
        idx = np.array([j for j in range(len(y)) if j != i])
        b = fit_beta(idx, **cfg)["beta"]
        pred[i] = Zt[i] @ b
    auc = roc_auc_score(y, pred)
    loo_rows.append({**cfg, "loo_auc": auc})
    print(f"  {cfg} -> LOO AUC {auc:.3f}   ({time.time()-t0:.0f}s)")
loo = pd.DataFrame(loo_rows).sort_values("loo_auc", ascending=False)
loo.to_csv(f"{OUT}/jse_loo_tuning.csv", index=False)
best = loo.iloc[0].to_dict()
cfg = {k: best[k] for k in ["tau0", "omega0", "gene_level", "module_level"]}
cfg["gene_level"] = bool(cfg["gene_level"]); cfg["module_level"] = bool(cfg["module_level"])
print("\nLOCKED configuration:", cfg, "JSE LOO AUC =", round(best["loo_auc"], 3))

# ---------------- 2. lock: refit on all 20 discovery patients ----------------
post = fit_beta(np.arange(len(y)), **cfg)
beta = post["beta"]
pd.DataFrame({"gene": genes, "beta": beta, "beta_sd": post["beta_sd"]}).sort_values(
    "beta", key=np.abs, ascending=False).to_csv(f"{OUT}/final_gene_weights.csv", index=False)
theta = pd.DataFrame({"module": mnames, "source": msrc, "theta": post["theta"],
                      "theta_sd": post["theta_sd"], "p_direction": post["theta_pdir"]})
theta.reindex(theta.theta.abs().sort_values(ascending=False).index).to_csv(f"{OUT}/final_module_weights.csv", index=False)
print("\ntop modules by |theta|:")
print(theta.reindex(theta.theta.abs().sort_values(ascending=False).index).head(10).to_string(index=False))

score = Zd @ beta
score_z = stats.zscore(score)
np.save(f"{OUT}/dk_final_score.npy", score)

# ---------------- 3. validation on DK ----------------
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
early, late = (T <= 24) & (E == 1), T > 60
bmask = early | late


def pack(T_, E_, sc, tag):
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hrs = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(hrs.confidence_intervals_.iloc[0].values)
    return dict(endpoint=tag, C=float(c), hr_per_sd=float(np.exp(cph.params_.iloc[0])),
                cox_p=float(cph.summary["p"].iloc[0]), logrank_p=float(lr.p_value),
                hr_median_split=float(np.exp(hrs.params_.iloc[0])),
                hr_ci=[float(ci[0]), float(ci[1])])


res = {"config": cfg, "jse_loo_auc": float(best["loo_auc"]),
       "primary_OS_5yr": pack(T5, E5, score, "OS 5-year (primary)"),
       "secondary_OS_full": pack(T, E, score, "OS full follow-up"),
       "binary_early_vs_long": dict(auc=float(roc_auc_score(early[bmask].astype(int), score[bmask])),
                                    n_early=int(early.sum()), n_long=int(late.sum()),
                                    mw_p=float(stats.mannwhitneyu(score[early], score[late]).pvalue))}
print("\nPRIMARY (5-year OS):", res["primary_OS_5yr"])
print("SECONDARY (full OS):", res["secondary_OS_full"])
print("BINARY early-death vs long-survivor:", res["binary_early_vs_long"])

# clinical adjustment
cl = pd.DataFrame({"time": T5, "event": E5, "score": score_z}, index=clin.index)
cl["stage_III_IV"] = (clin["STAGE"].astype(str) == "III-IV").astype(int)
cl["solid"] = (clin["FORM"].astype(str) == "Solid").astype(int)
u = clin["FORM"].astype(str).isin(["Solid", "Tubulocribriform"])
mm = CoxPHFitter().fit(cl.loc[u, ["time", "event", "score", "stage_III_IV", "solid"]], "time", "event")
res["multivariable_cox_5yr"] = {k: dict(HR=float(np.exp(mm.params_[k])), p=float(mm.summary.loc[k, "p"]))
                                for k in ["score", "stage_III_IV", "solid"]}
print("\nmultivariable (5-yr OS):\n", mm.summary[["exp(coef)", "p"]].to_string())

# nulls
c_obs = res["primary_OS_5yr"]["C"]
nz = np.count_nonzero(np.abs(beta) > np.percentile(np.abs(beta), 90))
rand = []
for b in range(2000):
    idx = rng.choice(len(genes), nz, replace=False)
    sgn = rng.choice([-1, 1], nz)
    rand.append(concordance_index_censored(E5.astype(bool), T5, Zd[:, idx] @ sgn)[0])
rand = np.array(rand)
perm = []
for b in range(2000):
    o = rng.permutation(len(T))
    perm.append(concordance_index_censored(E5[o].astype(bool), T5[o], score)[0])
perm = np.array(perm)
res["nulls"] = dict(random_signature_p=float((np.sum(rand >= c_obs) + 1) / 2001),
                    permutation_p=float((np.sum(np.array(perm) >= c_obs) + 1) / 2001),
                    random_mean=float(rand.mean()), perm_mean=float(perm.mean()))
print("\nnulls:", res["nulls"])

# ---------------- 4. figures ----------------
def km(sc, T_, E_, title, fname):
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hr = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event").params_.iloc[0])
    fig, ax = plt.subplots(figsize=(5.4, 4.3))
    for m_, lab, col in [(hi, "high risk", "#c0392b"), (~hi, "low risk", "#2471a3")]:
        KaplanMeierFitter().fit(T_[m_], E_[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, color=col)
    ax.set_xlabel("months"); ax.set_ylabel("overall survival"); ax.set_ylim(0, 1.02)
    ax.set_title(f"{title}\nHR = {hr:.2f}, log-rank p = {lr.p_value:.4f}", fontsize=10)
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/{fname}", dpi=150); plt.close(fig)


km(score, T5, E5, "DK external cohort — 5-year OS (primary)", "FINAL_km_5yr.png")
km(score, T, E, "DK external cohort — full follow-up", "FINAL_km_full.png")

q = np.quantile(score, [1/3, 2/3]); grp = np.digitize(score, q)
lr3 = multivariate_logrank_test(T5, grp, E5)
fig, ax = plt.subplots(figsize=(5.6, 4.3))
for g_, lab, col in [(2, "high", "#c0392b"), (1, "intermediate", "#7f8c8d"), (0, "low", "#2471a3")]:
    m_ = grp == g_
    KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, ci_show=False, color=col)
ax.set_xlabel("months"); ax.set_ylabel("overall survival"); ax.set_ylim(0, 1.02)
ax.set_title(f"DK external cohort — signature tertiles, 5-year OS\nlog-rank p = {lr3.p_value:.4f}", fontsize=10)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/FINAL_km_tertiles.png", dpi=150); plt.close(fig)
hi3, lo3 = grp == 2, grp == 0
lrx = logrank_test(T5[hi3], T5[lo3], E5[hi3], E5[lo3])
res["tertiles_5yr"] = dict(logrank_p_3grp=float(lr3.p_value), logrank_p_top_vs_bottom=float(lrx.p_value))
print("tertiles:", res["tertiles_5yr"])

# ---------------- 5. compact panel from the locked weights ----------------
w = pd.Series(beta, index=genes)
top = w.reindex(w.abs().sort_values(ascending=False).index)
for npan in [10, 20, 30, 50]:
    g = top.index[:npan]
    sc = (Zd[:, [genes.index(x) for x in g]] * np.sign(top[:npan].values)).mean(1)
    r = pack(T5, E5, sc, f"panel{npan}")
    res[f"panel{npan}"] = {**r, "genes": list(g)}
    print(f"panel {npan}: C={r['C']:.3f} logrank p={r['logrank_p']:.4g} HR={r['hr_median_split']:.2f}")
km((Zd[:, [genes.index(x) for x in top.index[:20]]] * np.sign(top[:20].values)).mean(1), T5, E5,
   "DK external cohort — 20-gene panel, 5-year OS", "FINAL_km_panel20.png")

json.dump(res, open(f"{OUT}/FINAL_validation.json", "w"), indent=2, default=float)
json.dump(dict(config=cfg, genes=genes, beta=beta.tolist()), open(f"{OUT}/FINAL_model.json", "w"))
print("\nDONE ->", OUT)

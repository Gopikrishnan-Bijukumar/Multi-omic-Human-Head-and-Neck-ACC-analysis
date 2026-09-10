"""Final locked signature, compact biomarker panel, and full validation package.

Outputs (work/outputs):
  final_signature.json        locked gene lists + scoring rule
  panel_genes.csv             compact panel with per-gene statistics in both cohorts
  km_*.png                    Kaplan-Meier splits
  validation_summary.csv      C-index, HR, logrank, permutation p, adjusted Cox, JSE AUC
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, json, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc
from sksurv.util import Surv
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD

rng = np.random.default_rng(6)
D = d3.load("all")
genes, T, E, y = D["genes"], D["T"], D["E"], D["y"]
Zd, Zt = D["Zd"].values, D["Zt"].values
clin = D["clin"]
M, mnames, msrc = MOD.build(genes, R)
SC = Zd @ M
K = M.shape[1]
gi = {g: i for i, g in enumerate(genes)}


def cox_z(X, T_, E_):
    o = np.argsort(T_); Xs = X[o]; Es = E_[o]
    Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
    n = Xs.shape[0]
    cs = np.cumsum(Xs[::-1], 0)[::-1]; cs2 = np.cumsum((Xs ** 2)[::-1], 0)[::-1]
    nr = np.arange(n, 0, -1)[:, None]
    mr = cs / nr; vr = cs2 / nr - mr ** 2
    m = Es == 1
    return (Xs[m] - mr[m]).sum(0) / np.sqrt(vr[m].sum(0) + 1e-12)


def topk2_score(tr, te, k=2):
    z = cox_z(SC[tr], T[tr], E[tr])
    top = np.argsort(-np.abs(z))[:k]
    zz = cox_z(SC[tr][:, top], T[tr], E[tr])
    mu, sd = SC[tr][:, top].mean(0), SC[tr][:, top].std(0) + 1e-9
    return (((SC[te][:, top] - mu) / sd) * np.sign(zz)).mean(1), top, np.sign(zz)


def km_plot(score, name, title, fname, labels=("high risk", "low risk")):
    hi = score > np.median(score)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    hr = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T, "event": E, "g": hi.astype(int)}),
                                  "time", "event").params_.iloc[0])
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for mask, lab, col in [(hi, labels[0], "#c0392b"), (~hi, labels[1], "#2471a3")]:
        KaplanMeierFitter().fit(T[mask], E[mask], label=f"{lab} (n={mask.sum()})").plot_survival_function(ax=ax, ci_show=True, color=col)
    ax.set_xlabel("months"); ax.set_ylabel("overall survival")
    ax.set_title(f"{title}\nHR={hr:.2f}, log-rank p={lr.p_value:.4f}", fontsize=10)
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/{fname}", dpi=150); plt.close(fig)
    return hr, lr.p_value


results = {}

# ============ 1. honest out-of-fold performance of the locked PROCEDURE ============
N_REP = 20
oof = np.zeros((N_REP, len(T)))
picks = []
for rep in range(N_REP):
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=200 + rep).split(Zd, E):
        sc, top, sg = topk2_score(tr, te)
        oof[rep, te] = stats.zscore(sc) if np.std(sc) > 1e-12 else 0
        picks.append(tuple(sorted(mnames[i] for i in top)))
oof_m = oof.mean(0)
c_oof = concordance_index_censored(E.astype(bool), T, oof_m)[0]
hr_oof, p_oof = km_plot(oof_m, "oof", "DK cohort — out-of-fold risk score\n(module selection inside every fold)", "km_oof_procedure.png")
results["oof_procedure"] = dict(C=c_oof, hr_split=hr_oof, logrank_p=p_oof,
                                c_per_rep_mean=float(np.mean([concordance_index_censored(E.astype(bool), T, oof[r])[0] for r in range(N_REP)])))
print("OOF procedure:", results["oof_procedure"])
print("most frequent module pairs:", pd.Series(picks).value_counts().head(5).to_dict())

# ============ 2. permutation null for the whole procedure ============
n_perm = 300
perm_c = []
for b in range(n_perm):
    o = rng.permutation(len(T))
    Tp, Ep = T[o], E[o]
    op = np.zeros(len(T))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=7).split(Zd, Ep):
        z = cox_z(SC[tr], Tp[tr], Ep[tr]); top = np.argsort(-np.abs(z))[:2]
        zz = cox_z(SC[tr][:, top], Tp[tr], Ep[tr])
        mu, sd = SC[tr][:, top].mean(0), SC[tr][:, top].std(0) + 1e-9
        op[te] = (((SC[te][:, top] - mu) / sd) * np.sign(zz)).mean(1)
    perm_c.append(concordance_index_censored(Ep.astype(bool), Tp, op)[0])
perm_c = np.array(perm_c)
p_perm = (np.sum(perm_c >= c_oof) + 1) / (n_perm + 1)
results["permutation"] = dict(n=n_perm, null_mean=float(perm_c.mean()), null_p95=float(np.percentile(perm_c, 95)), p=float(p_perm))
print("permutation null:", results["permutation"])

# ============ 3. locked signature on the full cohort ============
z_full = cox_z(SC, T, E)
top_full = np.argsort(-np.abs(z_full))[:2]
sel_modules = [mnames[i] for i in top_full]
sel_sign = np.sign(z_full[top_full])
print("\nselected modules (full DK):", list(zip(sel_modules, sel_sign)), "z =", z_full[top_full].round(2))
sig_score = (((SC[:, top_full] - SC[:, top_full].mean(0)) / (SC[:, top_full].std(0) + 1e-9)) * sel_sign).mean(1)
hr_lock, p_lock = km_plot(sig_score, "locked", f"DK cohort — locked signature ({' + '.join(sel_modules)})", "km_locked_signature.png")
c_lock = concordance_index_censored(E.astype(bool), T, sig_score)[0]
results["locked_full"] = dict(modules=sel_modules, signs=sel_sign.tolist(), C=c_lock, hr_split=hr_lock, logrank_p=p_lock)

# tertile split (top third vs bottom third) — the clinically usable presentation
q = np.quantile(sig_score, [1/3, 2/3])
grp = np.digitize(sig_score, q)
lr3 = multivariate_logrank_test(T, grp, E)
fig, ax = plt.subplots(figsize=(5.6, 4.4))
for g, lab, col in [(2, "high risk", "#c0392b"), (1, "intermediate", "#7f8c8d"), (0, "low risk", "#2471a3")]:
    m = grp == g
    KaplanMeierFitter().fit(T[m], E[m], label=f"{lab} (n={m.sum()})").plot_survival_function(ax=ax, ci_show=False, color=col)
ax.set_xlabel("months"); ax.set_ylabel("overall survival")
ax.set_title(f"DK cohort — signature tertiles\nlog-rank p={lr3.p_value:.4f}", fontsize=10)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/km_tertiles.png", dpi=150); plt.close(fig)
hi3, lo3 = grp == 2, grp == 0
lr_ext = logrank_test(T[hi3], T[lo3], E[hi3], E[lo3])
hr_ext = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": np.r_[T[hi3], T[lo3]], "event": np.r_[E[hi3], E[lo3]],
                                                "g": np.r_[np.ones(hi3.sum()), np.zeros(lo3.sum())]}), "time", "event").params_.iloc[0])
results["tertiles"] = dict(logrank_p_3grp=float(lr3.p_value), hr_top_vs_bottom=float(hr_ext), logrank_p_top_vs_bottom=float(lr_ext.p_value))
print("tertiles:", results["tertiles"])

# ============ 4. random-signature null (is it the genes, or any gene set?) ============
sizes = [int(np.sum(M[:, i] != 0)) for i in top_full]
rand_c = []
for b in range(2000):
    s = np.zeros(len(T))
    for sz, sg in zip(sizes, sel_sign):
        idx = rng.choice(len(genes), sz, replace=False)
        s += stats.zscore(Zd[:, idx].mean(1)) * sg
    rand_c.append(concordance_index_censored(E.astype(bool), T, s)[0])
rand_c = np.array(rand_c)
results["random_signature_null"] = dict(n=2000, p=float((np.sum(rand_c >= c_lock) + 1) / 2001),
                                        null_mean=float(rand_c.mean()))
print("random-signature null:", results["random_signature_null"])

# ============ 5. independence from clinical stage / histology ============
cl = pd.DataFrame({"time": T, "event": E, "score": stats.zscore(sig_score)}, index=clin.index)
cl["stage_III_IV"] = (clin["STAGE"].astype(str) == "III-IV").astype(int)
cl["solid"] = (clin["FORM"].astype(str) == "Solid").astype(int)
usable = clin["FORM"].astype(str).isin(["Solid", "Tubulocribriform"])
m_uni = CoxPHFitter().fit(cl[["time", "event", "score"]], "time", "event")
m_multi = CoxPHFitter().fit(cl.loc[usable, ["time", "event", "score", "stage_III_IV", "solid"]], "time", "event")
print("\nmultivariable Cox (score + stage + histology):\n", m_multi.summary[["coef", "exp(coef)", "p"]].to_string())
results["multivariable_cox"] = {k: dict(HR=float(np.exp(m_multi.params_[k])), p=float(m_multi.summary.loc[k, "p"]))
                                for k in ["score", "stage_III_IV", "solid"]}
results["univariable_cox"] = dict(HR=float(np.exp(m_uni.params_["score"])), p=float(m_uni.summary.loc["score", "p"]))

# ============ 6. compact gene panel ============
# leading-edge: genes of the selected modules ranked by their own DK association,
# with the JSE discovery direction reported alongside (never used to filter here).
mem = np.abs(M[:, top_full]).sum(1) > 0
cand_idx = np.where(mem)[0]
gz = cox_z(Zd[:, cand_idx], T, E)
tj = stats.ttest_ind(Zt[:, cand_idx][y == 1], Zt[:, cand_idx][y == 0], axis=0)
panel = pd.DataFrame({
    "gene": [genes[i] for i in cand_idx],
    "module": ["/".join([sel_modules[j] for j in range(len(top_full)) if M[i, top_full[j]] != 0]) for i in cand_idx],
    "dk_cox_z": gz,
    "dk_HR_per_SD": np.exp(gz / np.sqrt(E.sum())),
    "jse_t_PoorvsGood": tj.statistic,
    "jse_p": tj.pvalue,
}).sort_values("dk_cox_z", ascending=False)
panel["direction_agreement"] = np.sign(panel.dk_cox_z) == np.sign(panel.jse_t_PoorvsGood)
panel.to_csv(f"{OUT}/panel_genes.csv", index=False)
print(f"\npanel candidates n={len(panel)}; direction agreement with JSE = {panel.direction_agreement.mean():.2f}")
print(panel.head(20).to_string(index=False))

# compact 20-gene panel, honestly cross-validated (gene picking inside folds)
def panel_cv(npanel):
    o = np.zeros((10, len(T)))
    for rep in range(10):
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=300 + rep).split(Zd, E):
            zf = cox_z(SC[tr], T[tr], E[tr]); tp2 = np.argsort(-np.abs(zf))[:2]
            memf = np.where(np.abs(M[:, tp2]).sum(1) > 0)[0]
            gzf = cox_z(Zd[tr][:, memf], T[tr], E[tr])
            pick = memf[np.argsort(-np.abs(gzf))[:npanel]]
            sgn = np.sign(gzf[np.argsort(-np.abs(gzf))[:npanel]])
            mu, sd = Zd[tr][:, pick].mean(0), Zd[tr][:, pick].std(0) + 1e-9
            o[rep, te] = (((Zd[te][:, pick] - mu) / sd) * sgn).mean(1)
    om = o.mean(0)
    return concordance_index_censored(E.astype(bool), T, om)[0], om

for npan in [10, 15, 20, 30, 50]:
    c, om = panel_cv(npan)
    hi = om > np.median(om)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    print(f"CV-honest compact panel n={npan}: C={c:.3f} logrank p={lr.p_value:.4g}")
    results[f"panel_cv_{npan}"] = dict(C=float(c), logrank_p=float(lr.p_value))

top20 = panel.reindex(panel.dk_cox_z.abs().sort_values(ascending=False).index).head(20)
p20_score = ((Zd[:, [gi[g] for g in top20.gene]] - Zd[:, [gi[g] for g in top20.gene]].mean(0)) /
             (Zd[:, [gi[g] for g in top20.gene]].std(0) + 1e-9) * np.sign(top20.dk_cox_z.values)).mean(1)
hr20, p20 = km_plot(p20_score, "panel20", "DK cohort — 20-gene compact panel (full-cohort fit)", "km_panel20.png")
results["panel20_full"] = dict(genes=top20.gene.tolist(), hr_split=float(hr20), logrank_p=float(p20),
                               C=float(concordance_index_censored(E.astype(bool), T, p20_score)[0]))

# ============ 7. external check in the JSE discovery cohort ============
sct = Zt @ M[:, top_full]
sig_jse = (((sct - sct.mean(0)) / (sct.std(0) + 1e-9)) * sel_sign).mean(1)
auc_jse = roc_auc_score(y, sig_jse)
results["jse_external"] = dict(auc_poor_vs_good=float(auc_jse),
                               mannwhitney_p=float(stats.mannwhitneyu(sig_jse[y == 1], sig_jse[y == 0]).pvalue))
print("\nJSE cohort (independent of the DK fit): AUC Poor vs Good =", round(auc_jse, 3), results["jse_external"])

# ============ 8. time-dependent AUC ============
try:
    sv = Surv.from_arrays(E.astype(bool), T)
    times = np.array([24, 60, 120])
    times = times[times < T[E == 1].max()]
    auc_t, mean_auc = cumulative_dynamic_auc(sv, sv, sig_score, times)
    results["time_dependent_auc"] = dict(times=times.tolist(), auc=auc_t.tolist(), mean=float(mean_auc))
    print("time-dependent AUC:", dict(zip(times, auc_t.round(3))))
except Exception as ex:
    print("tdAUC failed", ex)

json.dump(results, open(f"{OUT}/validation_summary.json", "w"), indent=2, default=float)
json.dump(dict(modules=sel_modules, signs=sel_sign.tolist(),
               gene_sets={m: [genes[i] for i in np.where(M[:, mnames.index(m)] != 0)[0]] for m in sel_modules},
               scoring="per-gene rank-inverse-normal within cohort; mean per module; z-score module scores; signed mean",
               panel20=top20.gene.tolist()),
          open(f"{OUT}/final_signature.json", "w"), indent=2)
print("\nwrote", OUT)

"""COMPACT PANEL: how few genes can carry the signature?

Nothing here overwrites the 20-gene result. All new files go to work/outputs_compact/.

Protocol (panel size and membership are chosen WITHOUT looking at DK)
--------------------------------------------------------------------
A. Resampling inside the discovery cohort only:
     - leave-one-out (20 folds) x 3 seeds
     - repeated stratified 5-fold (5 repeats x 5 folds)
   In every fold the model is refitted on the training patients, the top-k genes by
   |beta| are taken, and the held-out patient(s) are scored with that fold's panel.
   -> honest discovery AUC as a function of panel size k, and gene selection frequency.
B. Pre-specified parsimony rule: choose the SMALLEST k whose discovery AUC is within
   0.02 of the best k. Scoring scheme (sign-mean vs beta-weighted) chosen the same way.
C. Lock: panel = top-k genes of the already-locked full-discovery beta (FINAL_model.json).
D. DK is touched once: C-index, HR, log-rank, tertiles, binary endpoint, permutation and
   random-gene-set nulls, clinical multivariable model.
E. Adoption rule: keep the compact panel only if it is not worse than the 20-gene panel
   on DK (C within 0.01 and log-rank p < 0.05). Otherwise the 20-gene panel stands.
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
from sklearn.model_selection import RepeatedStratifiedKFold
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_compact"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
import modules as MOD
from mlbayes import MultiLevelModel, to_t

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})


def save(fig, name):
    """Every figure in both formats; PNG at 600 dpi."""
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


rng = np.random.default_rng(11)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]

keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gi = {g: i for i, g in enumerate(genes)}
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
M, mnames, msrc = MOD.build(genes, R)
srcs = sorted(set(msrc)); src_ids = np.array([srcs.index(s) for s in msrc])
ev_tab = pd.read_csv(f"{R}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv").set_index("gene")
s_ev, d_ev = d3.evidence_prior(ev_tab.reindex(genes))

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
cfg = LOCK["config"]
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
assert not np.isnan(beta_lock).any()
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
print(f"genes={len(genes)} modules={M.shape[1]} JSE n={len(y)} DK n={len(T)} events5y={E5.sum()}")
print("locked config:", cfg)

KS = [3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30]
STEPS = 900


def fit_beta(idx, seed=0):
    m = MultiLevelModel(M, src_ids, len(srcs), s_ev, d_ev, tau0=cfg["tau0"], omega0=cfg["omega0"],
                        use_gene_level=bool(cfg["gene_level"]), use_module_level=bool(cfg["module_level"]),
                        use_dk=False, use_jse=True, seed=seed)
    m.fit(Xj=to_t(Zt[idx]), yj=to_t(y[idx]), steps=STEPS)
    return m.posterior(300)["beta"]


def panel_score(Z, idx_g, w, scheme):
    """Panel score for rows of Z using gene positions idx_g and their training weights w."""
    if scheme == "sign":
        return (Z[:, idx_g] * np.sign(w)).mean(1)
    return Z[:, idx_g] @ w


# ---------------- A. discovery-internal resampling ----------------
t0 = time.time()
folds = []
for seed in range(3):                                    # leave-one-out x 3 seeds
    for i in range(len(y)):
        folds.append((np.array([j for j in range(len(y)) if j != i]), np.array([i]), seed, "loo"))
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=7)
for f, (tr, te) in enumerate(rskf.split(np.zeros(len(y)), y)):
    folds.append((tr, te, 100 + f, "cv5"))
print(f"{len(folds)} training fits")

pred = {("loo", k, s): np.full((3, len(y)), np.nan) for k in KS for s in ("sign", "beta")}
pred_cv = {("cv5", k, s): [] for k in KS for s in ("sign", "beta")}
sel_count = np.zeros((len(KS), len(genes)))
sign_count = np.zeros(len(genes)); sign_tot = 0

for n, (tr, te, seed, kind) in enumerate(folds):
    b = fit_beta(tr, seed=seed)
    order = np.argsort(-np.abs(b))
    sign_count += np.sign(b) * (np.abs(b) > np.percentile(np.abs(b), 99)); sign_tot += 1
    for ki, k in enumerate(KS):
        idx_g = order[:k]
        sel_count[ki, idx_g] += 1
        for scheme in ("sign", "beta"):
            sc = panel_score(Zt[te], idx_g, b[idx_g], scheme)
            if kind == "loo":
                pred[("loo", k, scheme)][seed, te[0]] = sc[0]
            else:
                pred_cv[("cv5", k, scheme)].append((te, sc))
    if (n + 1) % 20 == 0:
        print(f"  {n+1}/{len(folds)} fits  ({time.time()-t0:.0f}s)")

rows = []
for k in KS:
    for scheme in ("sign", "beta"):
        aucs_loo = [roc_auc_score(y, pred[("loo", k, scheme)][s]) for s in range(3)]
        cv = pred_cv[("cv5", k, scheme)]
        aucs_cv = []
        for rep in range(5):
            p = np.full(len(y), np.nan)
            for te, sc in cv[rep * 5:(rep + 1) * 5]:
                p[te] = sc
            aucs_cv.append(roc_auc_score(y, p))
        rows.append(dict(k=k, scheme=scheme, loo_auc=np.mean(aucs_loo), loo_sd=np.std(aucs_loo),
                         cv5_auc=np.mean(aucs_cv), cv5_sd=np.std(aucs_cv),
                         disc_auc=0.5 * np.mean(aucs_loo) + 0.5 * np.mean(aucs_cv)))
disc = pd.DataFrame(rows)
disc.to_csv(f"{OUT}/discovery_panel_size_scan.csv", index=False)
print("\ndiscovery-internal AUC by panel size (DK never seen):")
print(disc.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

# ---------------- B. pre-specified parsimony rule ----------------
best = disc["disc_auc"].max()
elig = disc[disc["disc_auc"] >= best - 0.02].sort_values(["k", "disc_auc"], ascending=[True, False])
K_STAR = int(elig.iloc[0]["k"]); SCHEME = elig.iloc[0]["scheme"]
print(f"\nbest discovery AUC = {best:.3f}; smallest k within 0.02 -> k={K_STAR}, scheme={SCHEME}")

# ---------------- C. lock the compact panel ----------------
order_lock = np.argsort(-np.abs(beta_lock))
panel_idx = order_lock[:K_STAR]
panel = [genes[i] for i in panel_idx]
freq = sel_count[KS.index(K_STAR)] / len(folds)
gene2mod = {}
for j, mn in enumerate(mnames):
    for i in np.where(M[:, j] != 0)[0]:
        gene2mod.setdefault(genes[i], []).append(mn)
ptab = pd.DataFrame({
    "gene": panel,
    "beta_locked": beta_lock[panel_idx],
    "direction": np.where(beta_lock[panel_idx] > 0, "poor-prognosis", "good-prognosis"),
    "selection_freq": freq[panel_idx],
    "evidence_direction": ev_tab.reindex(panel)["direction_consensus"].values,
    "n_evidence_sources": ev_tab.reindex(panel)["n_sources"].values if "n_sources" in ev_tab.columns else np.nan,
    "modules": [";".join(gene2mod.get(g, [])[:6]) for g in panel]})
ptab.to_csv(f"{OUT}/compact_panel_genes.csv", index=False)
print("\nlocked compact panel:")
print(ptab[["gene", "beta_locked", "direction", "selection_freq"]].to_string(index=False))

sc_compact = panel_score(Zd, panel_idx, beta_lock[panel_idx], SCHEME)
np.save(f"{OUT}/dk_compact_score.npy", sc_compact)
idx20 = order_lock[:20]
sc20 = panel_score(Zd, idx20, beta_lock[idx20], "sign")
sc_full = Zd @ beta_lock


# ---------------- D. DK validation ----------------
def pack(sc, T_, E_, tag):
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


res = {"locked_config": cfg, "chosen_k": K_STAR, "scheme": SCHEME,
       "discovery_best_auc": float(best), "panel": panel}
res["compact_5yr"] = pack(sc_compact, T5, E5, f"{K_STAR}-gene, 5-year OS (primary)")
res["compact_full"] = pack(sc_compact, T, E, f"{K_STAR}-gene, full follow-up")
res["panel20_5yr"] = pack(sc20, T5, E5, "20-gene, 5-year OS")
res["full_signature_5yr"] = pack(sc_full, T5, E5, "full signature, 5-year OS")
for k, v in [("compact_5yr", res["compact_5yr"]), ("compact_full", res["compact_full"]),
             ("panel20_5yr", res["panel20_5yr"]), ("full_signature_5yr", res["full_signature_5yr"])]:
    print(f"\n{v['endpoint']}: C={v['C']:.3f} HR/SD={v['hr_per_sd']:.2f} (p={v['cox_p']:.4f}) "
          f"split HR={v['hr_median_split']:.2f} [{v['hr_ci'][0]:.2f}-{v['hr_ci'][1]:.2f}] logrank p={v['logrank_p']:.4f}")

early, late = (T <= 24) & (E == 1), T > 60
bmask = early | late
res["binary_early_vs_long"] = dict(auc=float(roc_auc_score(early[bmask].astype(int), sc_compact[bmask])),
                                   mw_p=float(stats.mannwhitneyu(sc_compact[early], sc_compact[late]).pvalue),
                                   n_early=int(early.sum()), n_long=int(late.sum()))
print("binary early-death vs long-survivor:", res["binary_early_vs_long"])

q = np.quantile(sc_compact, [1 / 3, 2 / 3]); grp = np.digitize(sc_compact, q)
res["tertiles_5yr"] = dict(logrank_p_3grp=float(multivariate_logrank_test(T5, grp, E5).p_value),
                           logrank_p_top_vs_bottom=float(logrank_test(T5[grp == 2], T5[grp == 0],
                                                                      E5[grp == 2], E5[grp == 0]).p_value))
print("tertiles:", res["tertiles_5yr"])

cl = pd.DataFrame({"time": T5, "event": E5, "score": stats.zscore(sc_compact)}, index=clin.index)
cl["stage_III_IV"] = (clin["STAGE"].astype(str) == "III-IV").astype(int)
cl["solid"] = (clin["FORM"].astype(str) == "Solid").astype(int)
u = clin["FORM"].astype(str).isin(["Solid", "Tubulocribriform"])
mm = CoxPHFitter().fit(cl.loc[u, ["time", "event", "score", "stage_III_IV", "solid"]], "time", "event")
res["multivariable_cox_5yr"] = {k: dict(HR=float(np.exp(mm.params_[k])), p=float(mm.summary.loc[k, "p"]))
                                for k in ["score", "stage_III_IV", "solid"]}
print("\nmultivariable (5-yr OS):\n", mm.summary[["exp(coef)", "p"]].to_string())

# nulls for the compact panel
c_obs = res["compact_5yr"]["C"]
randc = np.array([concordance_index_censored(
    E5.astype(bool), T5, panel_score(Zd, rng.choice(len(genes), K_STAR, replace=False),
                                     rng.choice([-1.0, 1.0], K_STAR), "sign"))[0] for _ in range(2000)])
permc = np.array([(lambda o: concordance_index_censored(E5[o].astype(bool), T5[o], sc_compact)[0])(
    rng.permutation(len(T))) for _ in range(2000)])
res["nulls"] = dict(random_signature_p=float((np.sum(randc >= c_obs) + 1) / 2001),
                    permutation_p=float((np.sum(permc >= c_obs) + 1) / 2001),
                    random_mean=float(randc.mean()), perm_mean=float(permc.mean()))
print("nulls:", res["nulls"])

# leave-one-gene-out robustness of the compact panel on DK
logo = []
for j, g in enumerate(panel):
    keep_i = [panel_idx[m] for m in range(K_STAR) if m != j]
    s = panel_score(Zd, np.array(keep_i), beta_lock[keep_i], SCHEME)
    r = pack(s, T5, E5, f"drop {g}")
    logo.append(dict(dropped=g, C=r["C"], logrank_p=r["logrank_p"], hr=r["hr_median_split"]))
logo = pd.DataFrame(logo).sort_values("C")
logo.to_csv(f"{OUT}/leave_one_gene_out.csv", index=False)
res["leave_one_gene_out"] = dict(C_min=float(logo.C.min()), C_max=float(logo.C.max()),
                                 p_max=float(logo.logrank_p.max()))
print("\nleave-one-gene-out on DK:\n", logo.to_string(index=False))

# DK curve across all k (transparency only - not used for selection)
dk_curve = []
for k in KS:
    idxk = order_lock[:k]
    for scheme in ("sign", "beta"):
        r = pack(panel_score(Zd, idxk, beta_lock[idxk], scheme), T5, E5, f"k={k}/{scheme}")
        dk_curve.append(dict(k=k, scheme=scheme, C=r["C"], logrank_p=r["logrank_p"], hr=r["hr_median_split"]))
dk_curve = pd.DataFrame(dk_curve)
dk_curve.to_csv(f"{OUT}/dk_panel_size_curve.csv", index=False)

# ---------------- E. adoption rule ----------------
adopt = (res["compact_5yr"]["C"] >= res["panel20_5yr"]["C"] - 0.01) and (res["compact_5yr"]["logrank_p"] < 0.05)
res["adopt_compact"] = bool(adopt)
print(f"\nADOPT COMPACT PANEL: {adopt}")

# ---------------- figures ----------------
def km(sc, T_, E_, title, fname, ylab="overall survival"):
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    hr = np.exp(cf.params_.iloc[0]); ci = np.exp(cf.confidence_intervals_.iloc[0].values)
    fig, ax = plt.subplots(figsize=(4.6, 3.9))
    for m_, lab, col in [(hi, "high risk", "#b2182b"), (~hi, "low risk", "#2166ac")]:
        KaplanMeierFitter().fit(T_[m_], E_[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, color=col, ci_alpha=0.12)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nHR = {hr:.2f} [{ci[0]:.2f}–{ci[1]:.2f}], log-rank p = {lr.p_value:.3f}", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    save(fig, fname)


km(sc_compact, T5, E5, f"DK external cohort — {K_STAR}-gene panel, 5-year OS", "COMPACT_km_5yr")
km(sc_compact, T, E, f"DK external cohort — {K_STAR}-gene panel, full follow-up", "COMPACT_km_full")

fig, ax = plt.subplots(figsize=(5.2, 3.6))
for scheme, col in [("sign", "#b2182b"), ("beta", "#2166ac")]:
    d_ = disc[disc.scheme == scheme]
    ax.errorbar(d_.k, d_.disc_auc, yerr=d_.loo_sd, marker="o", ms=4, color=col, label=f"{scheme}-weighted", capsize=2)
ax.axvline(K_STAR, color="#555", ls="--", lw=1)
ax.annotate(f"chosen k = {K_STAR}", (K_STAR, ax.get_ylim()[0]), xytext=(4, 6),
            textcoords="offset points", fontsize=8, color="#555")
ax.set_xlabel("panel size (genes)"); ax.set_ylabel("discovery-internal AUC")
ax.set_title("Panel size chosen inside the discovery cohort\n(LOO + repeated 5-fold; DK not used)", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "COMPACT_panel_size_selection")

fig, ax = plt.subplots(figsize=(5.0, 0.32 * K_STAR + 1.4))
o = np.argsort(ptab.beta_locked.values)
cols = ["#b2182b" if b > 0 else "#2166ac" for b in ptab.beta_locked.values[o]]
ax.barh(range(K_STAR), ptab.beta_locked.values[o], color=cols)
ax.set_yticks(range(K_STAR)); ax.set_yticklabels(ptab.gene.values[o])
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("locked posterior mean weight (log-hazard per SD)")
ax.set_title(f"{K_STAR}-gene panel — red = poor prognosis, blue = protective", fontsize=9)
save(fig, "COMPACT_gene_weights")

fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
for ax, (vals, obs, lab) in zip(axes, [(randc, c_obs, f"random {K_STAR}-gene signatures"),
                                       (permc, c_obs, "permuted DK survival")]):
    ax.hist(vals, bins=40, color="#bbbbbb", edgecolor="none")
    ax.axvline(obs, color="#b2182b", lw=1.6)
    ax.set_xlabel("C-index (DK, 5-year OS)"); ax.set_ylabel("count"); ax.set_title(lab, fontsize=9)
axes[0].annotate(f"observed {c_obs:.3f}\np = {res['nulls']['random_signature_p']:.3f}",
                 (c_obs, axes[0].get_ylim()[1] * 0.8), xytext=(-70, 0), textcoords="offset points", fontsize=8, color="#b2182b")
axes[1].annotate(f"p = {res['nulls']['permutation_p']:.3f}", (c_obs, axes[1].get_ylim()[1] * 0.8),
                 xytext=(-58, 0), textcoords="offset points", fontsize=8, color="#b2182b")
save(fig, "COMPACT_nulls")

fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6), sharey=True)
for ax, (sc, name) in zip(axes, [(sc_full, "full signature"), (sc20, "20-gene panel"),
                                 (sc_compact, f"{K_STAR}-gene panel")]):
    hi = sc > np.median(sc)
    lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
    hr = np.exp(CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": hi.astype(int)}), "time", "event").params_.iloc[0])
    c = concordance_index_censored(E5.astype(bool), T5, sc)[0]
    for m_, lab, col in [(hi, "high risk", "#b2182b"), (~hi, "low risk", "#2166ac")]:
        KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(ax=ax, color=col, ci_alpha=0.1)
    ax.set_ylim(0, 1.02); ax.set_xlabel("months")
    ax.set_title(f"{name}\nC={c:.3f}  HR={hr:.2f}  p={lr.p_value:.3f}", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
axes[0].set_ylabel("overall survival (5-year)")
save(fig, "COMPACT_vs_full")

json.dump(res, open(f"{OUT}/COMPACT_validation.json", "w"), indent=2, default=float)
json.dump(dict(k=K_STAR, scheme=SCHEME, genes=panel,
               weights=dict(zip(panel, beta_lock[panel_idx].tolist()))),
          open(f"{OUT}/COMPACT_model.json", "w"), indent=2)
print("\nDONE ->", OUT)

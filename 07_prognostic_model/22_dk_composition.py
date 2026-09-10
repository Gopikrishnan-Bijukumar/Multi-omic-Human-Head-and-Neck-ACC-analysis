"""STAGE 4c: deconvolve DK and test the composition axis against survival.

Pre-specified in work/PRESPEC_stage4.md, written before this script was run and before any DK
survival was touched. The composition score's direction comes entirely from the TRUE single-cell
composition of the JSE patients:

    composition_score = z(Epithelial_Tumor) + z(Dividing) - z(Myoepithelial_Tumor)

Only rank/within-cohort-standardised quantities are used, because the sanity gate validated
cross-patient ordering (G1/G2 passed) but not absolute proportions (G3 failed).

The headline number is delta-C for adding composition to the 4-gene panel, with its CI.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import os, json
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import nnls
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")
import sys, importlib.util

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_deconv"
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY, GREEN, PURPLE = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837", "#762a83"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


META = json.load(open(f"{OUT}/reference_meta.json"))
COARSE = META["coarse_map"]
GATE = json.load(open(f"{OUT}/GATE_results.json"))
METHOD = GATE["best_method"]
print(f"deconvolution method (fixed by the gate): {METHOD}")

pb = pd.read_parquet(f"{OUT}/pseudobulk_ct_by_patient.parquet")
tr_log = pd.read_csv(f"{R}/outputs/03_harmonized/train_norm_log2.csv", index_col=0)
dk_log = pd.read_csv(f"{R}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
shared = sorted(set(pb.columns) & set(tr_log.columns) & set(dk_log.columns))
pb = pb[shared]
pb.index = pd.MultiIndex.from_arrays(
    [[COARSE.get(c, "Other") for c in pb.index.get_level_values("cell_type")],
     pb.index.get_level_values("patient")], names=["cell_type", "patient"])
pb = pb.groupby(level=["cell_type", "patient"]).sum()


def build_signature(n_markers=60):
    prof = pb.groupby(level="cell_type").sum()
    prof = prof.div(prof.sum(axis=1), axis=0) * 1e6
    lv = np.log2(prof + 1.0)
    markers = []
    for ct in prof.index:
        spec = lv.loc[ct] - lv.drop(index=ct).max(axis=0)
        markers += list(spec[prof.loc[ct] > 20].sort_values(ascending=False).index[:n_markers])
    markers = sorted(set(markers))
    return prof.T.loc[markers]


def dwls(bulk_vec, S):
    b = bulk_vec.values.astype(float); A = S.values.astype(float)
    w = nnls(A, b)[0]
    for _ in range(12):
        fit = A @ w
        wt = 1.0 / np.maximum(fit, np.percentile(fit[fit > 0], 5)) ** 2
        sw = np.sqrt(wt)
        w2 = nnls(A * sw[:, None], b * sw)[0]
        if np.max(np.abs(w2 - w)) < 1e-8:
            w = w2; break
        w = w2
    w = np.clip(w, 0, None)
    return pd.Series(w / w.sum() if w.sum() > 0 else np.full(len(w), 1 / len(w)), index=S.columns)


S = build_signature()
print(f"signature: {S.shape[0]} marker genes x {S.shape[1]} cell types")

# ---- survival + expression score (read-only) ----
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
gsub = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[gsub]).values
LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(gsub).values
idx4 = np.argsort(-np.abs(beta_lock))[:4]
expr4 = (Zd[:, idx4] * np.sign(beta_lock[idx4])).mean(1)
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
dk_ids = list(dk_raw.index)
n = len(T)
print(f"DK n={n}, events(5y)={E5.sum()}")

# ---- deconvolve DK (and JSE, as a labelled positive control) ----
dk_lin = (2.0 ** dk_log.loc[dk_ids, shared]) - 1.0; dk_lin[dk_lin < 0] = 0
tr_lin = (2.0 ** tr_log[shared]) - 1.0; tr_lin[tr_lin < 0] = 0
dk_comp = pd.DataFrame({p: dwls(dk_lin.loc[p, S.index], S) for p in dk_ids}).T
jse_comp = pd.DataFrame({p: dwls(tr_lin.loc[p, S.index], S) for p in tr_lin.index}).T
dk_comp.to_csv(f"{OUT}/dk_estimated_composition.csv")
jse_comp.to_csv(f"{OUT}/jse_estimated_composition_full.csv")

TRIAD = {"Epithelial_Tumor": +1.0, "Dividing": +1.0, "Myoepithelial_Tumor": -1.0}


def comp_score(df):
    z = (df - df.mean(0)) / (df.std(0) + 1e-12)
    return sum(w * z[c] for c, w in TRIAD.items()).values


comp_dk = comp_score(dk_comp)
comp_jse = comp_score(jse_comp)

res = {"method": METHOD, "gate_deviation": "G3 failed (0.685 vs 0.70); rank-based analysis only",
       "triad": TRIAD, "n_markers": int(S.shape[0])}

# positive control: does the score separate JSE Poor/Good? (direction came from JSE sc, so
# this is a consistency check, not independent evidence - labelled as such)
auc_jse = float(stats.mannwhitneyu(comp_jse[y == 1], comp_jse[y == 0]).pvalue)
from sklearn.metrics import roc_auc_score
res["jse_positive_control"] = dict(auc=float(roc_auc_score(y, comp_jse)), mw_p=auc_jse,
                                   note="direction derived from JSE single-cell; consistency check only")
print(f"\nJSE positive control (not independent): AUC={roc_auc_score(y, comp_jse):.3f}, p={auc_jse:.4f}")

# =====================================================================
# PRIMARY: Cox on continuous composition score, DK 5-year OS
# =====================================================================
print("\n" + "=" * 68 + "\nPRIMARY — composition score, DK 5-year OS\n" + "=" * 68)
cph = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "s": stats.zscore(comp_dk)}),
                        "time", "event")
c_comp = concordance_index_censored(E5.astype(bool), T5, comp_dk)[0]
hi = comp_dk > np.median(comp_dk)
lr = logrank_test(T5[hi], T5[~hi], E5[hi], E5[~hi])
res["primary"] = dict(hr_per_sd=float(np.exp(cph.params_.iloc[0])),
                      cox_p=float(cph.summary["p"].iloc[0]), C=float(c_comp),
                      logrank_p_median=float(lr.p_value),
                      significant=bool(cph.summary["p"].iloc[0] < 0.05))
print(f"HR/SD = {np.exp(cph.params_.iloc[0]):.3f}  p = {cph.summary['p'].iloc[0]:.4f}  "
      f"C = {c_comp:.3f}  median-split log-rank p = {lr.p_value:.4f}")

# =====================================================================
# SECONDARY: each compartment, FDR-corrected
# =====================================================================
rows = []
for ct in dk_comp.columns:
    v = dk_comp[ct].values
    if v.std() < 1e-9:
        continue
    cf = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "s": stats.zscore(v)}),
                           "time", "event")
    rows.append(dict(cell_type=ct, hr_per_sd=float(np.exp(cf.params_.iloc[0])),
                     p=float(cf.summary["p"].iloc[0]),
                     C=float(concordance_index_censored(E5.astype(bool), T5, v)[0])))
sec = pd.DataFrame(rows)
sec["fdr_q"] = multipletests(sec.p.values, method="fdr_bh")[1]
sec = sec.sort_values("p")
sec.to_csv(f"{OUT}/dk_percompartment_cox.csv", index=False)
res["secondary_percompartment"] = sec.to_dict("records")
print("\nSECONDARY — each compartment (FDR-corrected):")
print(sec.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# =====================================================================
# ORTHOGONALITY — the question Stage 4 exists to answer
# =====================================================================
rho, rho_p = stats.spearmanr(comp_dk, expr4)
res["orthogonality"] = dict(spearman_rho=float(rho), p=float(rho_p),
                            interpretation="low |rho| => genuinely different measurement axis")
print(f"\nORTHOGONALITY vs the 4-gene expression panel: Spearman rho = {rho:+.3f} (p = {rho_p:.4f})")

# =====================================================================
# COMBINATION — delta-C is the headline
# =====================================================================
bi = pd.DataFrame({"time": T5, "event": E5,
                   "expr4": stats.zscore(expr4), "comp": stats.zscore(comp_dk)})
bcox = CoxPHFitter().fit(bi, "time", "event")
res["bivariate_cox"] = {k: dict(HR=float(np.exp(bcox.params_[k])), p=float(bcox.summary.loc[k, "p"]))
                        for k in ["expr4", "comp"]}
print("\nBIVARIATE Cox:\n", bcox.summary[["exp(coef)", "p"]].to_string())

sgn = 1.0 if np.exp(bcox.params_["comp"]) >= 1 else -1.0
combined = stats.zscore(expr4) + sgn * stats.zscore(comp_dk)
c_expr = concordance_index_censored(E5.astype(bool), T5, expr4)[0]
c_comb = concordance_index_censored(E5.astype(bool), T5, combined)[0]
cph_c = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "s": stats.zscore(combined)}),
                          "time", "event")

rng = np.random.default_rng(41)
dC, dC_vs_comp = [], []
for _ in range(2000):
    b = rng.integers(0, n, n)
    if E5[b].sum() < 5:
        continue
    cc = concordance_index_censored(E5[b].astype(bool), T5[b], combined[b])[0]
    ce = concordance_index_censored(E5[b].astype(bool), T5[b], expr4[b])[0]
    co = concordance_index_censored(E5[b].astype(bool), T5[b], comp_dk[b])[0]
    dC.append(cc - ce); dC_vs_comp.append(cc - co)
res["combination"] = dict(
    C_expr4=float(c_expr), C_composition=float(c_comp), C_combined=float(c_comb),
    combined_hr_per_sd=float(np.exp(cph_c.params_.iloc[0])),
    combined_cox_p=float(cph_c.summary["p"].iloc[0]),
    deltaC_vs_expr4=dict(mean=float(np.mean(dC)), lo=float(np.percentile(dC, 2.5)),
                         hi=float(np.percentile(dC, 97.5))),
    deltaC_vs_composition=dict(mean=float(np.mean(dC_vs_comp)),
                               lo=float(np.percentile(dC_vs_comp, 2.5)),
                               hi=float(np.percentile(dC_vs_comp, 97.5))))
print(f"\nCOMBINATION:  C(expr4)={c_expr:.3f}  C(composition)={c_comp:.3f}  C(combined)={c_comb:.3f}")
print(f"  delta-C vs expression alone : {np.mean(dC):+.3f} "
      f"[{np.percentile(dC,2.5):+.3f}, {np.percentile(dC,97.5):+.3f}]   <-- HEADLINE")

# =====================================================================
# NULLS
# =====================================================================
perm = np.array([concordance_index_censored(E5[o].astype(bool), T5[o], comp_dk)[0]
                 for o in (rng.permutation(n) for _ in range(2000))])
zc = (dk_comp - dk_comp.mean(0)) / (dk_comp.std(0) + 1e-12)
randc = []
for _ in range(2000):
    cols = rng.choice(dk_comp.shape[1], 3, replace=False)
    sg = rng.choice([-1.0, 1.0], 3)
    randc.append(concordance_index_censored(E5.astype(bool), T5, (zc.values[:, cols] * sg).sum(1))[0])
randc = np.array(randc)
res["nulls"] = dict(permutation_p=float((np.sum(perm >= c_comp) + 1) / 2001),
                    random_triad_p=float((np.sum(randc >= c_comp) + 1) / 2001),
                    perm_mean=float(perm.mean()), rand_mean=float(randc.mean()))
print(f"\nNULLS: permutation p = {res['nulls']['permutation_p']:.4f}   "
      f"random 3-compartment p = {res['nulls']['random_triad_p']:.4f}")

json.dump(res, open(f"{OUT}/DK_COMPOSITION_results.json", "w"), indent=2, default=float)

# =====================================================================
# figures
# =====================================================================
fig, ax = plt.subplots(figsize=(5.2, 4.0))
early, late = (T <= 24) & (E == 1), T > 60
ax.scatter(stats.zscore(expr4)[~(early | late)], stats.zscore(comp_dk)[~(early | late)],
           c=GREY, s=30, label="intermediate", zorder=2)
ax.scatter(stats.zscore(expr4)[late], stats.zscore(comp_dk)[late], c=BLUE, s=42, label="alive >60 mo", zorder=3)
ax.scatter(stats.zscore(expr4)[early], stats.zscore(comp_dk)[early], c=RED, s=50, label="died <24 mo", zorder=4)
ax.axhline(0, color="k", lw=.6); ax.axvline(0, color="k", lw=.6)
ax.set_xlabel("4-gene expression score (z)"); ax.set_ylabel("cell-composition score (z)")
ax.set_title(f"Are expression and composition different axes?\nSpearman ρ = {rho:+.2f} (p = {rho_p:.3f})",
             fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "DK_composition_vs_expression")

fig, ax = plt.subplots(figsize=(4.7, 3.9))
for m_, lab, col in [(hi, "high composition score", RED), (~hi, "low", BLUE)]:
    KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{lab} (n={m_.sum()})").plot_survival_function(
        ax=ax, color=col, ci_alpha=.12)
ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("overall survival (5-year)")
ax.set_title(f"DK — cell-composition score\nHR/SD = {np.exp(cph.params_.iloc[0]):.2f}, "
             f"Cox p = {cph.summary['p'].iloc[0]:.3f}, log-rank p = {lr.p_value:.3f}", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "DK_composition_km")

f_ = sec.sort_values("hr_per_sd")
fig, ax = plt.subplots(figsize=(5.6, 0.35 * len(f_) + 1.4))
cols = [RED if p < 0.05 else GREY for p in f_.p.values]
ax.barh(range(len(f_)), np.log2(f_.hr_per_sd.values), color=cols)
ax.set_yticks(range(len(f_))); ax.set_yticklabels(f_.cell_type.values, fontsize=8)
ax.axvline(0, color="k", lw=.9)
ax.set_xlabel("log2 hazard ratio per SD of estimated fraction")
ax.set_title("DK survival by cell-type compartment\nred = nominal p < 0.05 (see FDR column)", fontsize=9)
save(fig, "DK_compartment_forest")

fig, ax = plt.subplots(figsize=(5.4, 3.5))
ax.bar(["expression\n(4 genes)", "composition\n(3 compartments)", "combined"],
       [c_expr - .5, c_comp - .5, c_comb - .5], bottom=.5,
       color=[RED, PURPLE, GREEN])
ax.axhline(.5, color="k", lw=.9)
for i, v in enumerate([c_expr, c_comp, c_comb]):
    ax.text(i, v + .004, f"{v:.3f}", ha="center", fontsize=9)
ax.set_ylim(.5, max(.75, c_comb + .05)); ax.set_ylabel("C-index (DK, 5-year OS)")
ax.set_title(f"Does composition add to expression?\nΔC = {np.mean(dC):+.3f} "
             f"[{np.percentile(dC,2.5):+.3f}, {np.percentile(dC,97.5):+.3f}]", fontsize=9)
save(fig, "DK_composition_deltaC")

print(f"\nDONE -> {OUT}")

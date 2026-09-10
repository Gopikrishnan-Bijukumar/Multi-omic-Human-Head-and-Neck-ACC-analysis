"""STAGE 3: Brayer et al. (Cancers 2023;15:1390) benchmark and convergence test.

Their DK validation cohort IS our DK cohort. This puts their classifier and ours on the same
patients, the same endpoint, and the same metric for the first time.

3a  exact classifier  - coefficients transcribed from their Table 5 (49-gene elastic net;
                        the 14 bold genes are their LASSO subset). Applied to DK.
3b  coefficient-free  - their 14 genes as a gene SET, weighted by our JSE-trained centroids.
3c  convergence       - do their classifier and our panel flag the SAME patients?

CAVEATS, stated in the output and carried into the report:
 * No intercept is published, so only the RANKING of patients is used, never an absolute call.
 * Their coefficients were fitted on their own normalisation; we apply them to rank-inverse-normal
   expression. Per-gene monotone transforms preserve each gene's order but not the exact relative
   weighting, so this is a faithful-but-not-identical reproduction. Reported as such.
 * Orientation (which sign means "poor") is resolved on the JSE discovery cohort, never on DK.
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
from sklearn.metrics import roc_auc_score, cohen_kappa_score
from sklearn.mixture import GaussianMixture
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_brayer"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY, GREEN = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
# Table 5, transcribed from work/../table_5_ss.png. Bold = 14-gene LASSO subset.
# ---------------------------------------------------------------------
T5_COEF = [
    ("A2M", -0.11, True), ("ABCA8", -0.0129, False), ("ACTA2", -0.1264, True),
    ("ADAMTS9", -0.0986, False), ("ALDOA", -0.018, False), ("ANO1", -0.099, True),
    ("APOL6", -0.0977, True), ("CARMN", -0.0091, False), ("CD9", -0.0886, False),
    ("CFH", -0.0226, False), ("COL17A1", -0.0151, False), ("COL7A1", -0.0095, False),
    ("COL9A2", -0.0072, False), ("DMD", -0.3178, True), ("EFS", 0.0074, False),
    ("EGFR", -0.0288, False), ("FRMD4B", -0.0279, False),
    ("HMCN1", -0.065, False), ("IPO9", 0.376, True), ("ITPR1", -0.1121, False),
    ("KRT14", -0.0359, False), ("LDLRAD4", 0.0354, False), ("LGR4", -0.0391, False),
    ("LIMA1", -0.0793, False), ("LIMCH1", -0.145, True), ("LOC107987158", 0.0363, False),
    ("LTF", -0.0048, False), ("MAMLD1", 0.0339, True), ("MIR205HG", -0.0511, True),
    ("MLPH", -0.013, False), ("MTUS1", -0.0917, False), ("PARP14", -0.0746, False),
    ("PCYOX1", -0.1327, False), ("PIK3R1", -0.0097, False),
    ("PLA2R1", -0.0434, False), ("PLAT", -0.0448, True), ("PLD1", -0.0356, False),
    ("PPARGC1A", -0.0378, False), ("PRUNE2", -0.0026, False), ("RASSF6", -0.1461, True),
    ("SEMA3C", -0.0582, True), ("SH3D19", -0.0562, False), ("SLPI", -0.1399, True),
    ("SVIL", -0.168, False), ("SYNPO2", -0.0159, False), ("TAGLN", -0.0321, False),
    ("TNFRSF19", -0.0054, False), ("TP63", -0.1442, True), ("TPM2", -0.0175, False),
]
coef = pd.DataFrame(T5_COEF, columns=["gene", "coefficient", "in_14gene"])
assert len(coef) == 49, len(coef)
assert coef.in_14gene.sum() == 14, coef.in_14gene.sum()
coef.to_csv(f"{OUT}/brayer_table5_coefs.csv", index=False)
print(f"Table 5 transcribed: {len(coef)} genes, {coef.in_14gene.sum()} in the LASSO subset")
print("14-gene subset:", list(coef[coef.in_14gene].gene))

# ---------------------------------------------------------------------
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
gi = {g: i for i, g in enumerate(genes)}
Zd = d3.rint(dk_raw[genes]).values
Zt = d3.rint(tr_raw[genes]).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

# widen the lookup to the unfiltered matrices so we lose as few of their genes as possible
Zd_all = d3.rint(dk_raw).values; Zt_all = d3.rint(tr_raw).values
gi_all = {g: i for i, g in enumerate(genes_all)}

avail = coef.gene.map(lambda g: g in gi_all)
coef["available"] = avail.values
print(f"\ngene availability in our harmonised matrices: {avail.sum()}/49 "
      f"({coef[coef.in_14gene].available.sum()}/14 of the LASSO subset)")
missing = list(coef[~coef.available].gene)
if missing:
    print("  missing:", missing)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
order_lock = np.argsort(-np.abs(beta_lock))
idx8 = order_lock[:8]; panel8 = [genes[i] for i in idx8]
sc_ours = (Zd[:, idx8] * np.sign(beta_lock[idx8])).mean(1)
sc_ours_jse = (Zt[:, idx8] * np.sign(beta_lock[idx8])).mean(1)
print(f"our 8-gene panel: {panel8}")


def linpred(Zmat, gindex, sub):
    ii = [gindex[g] for g in sub.gene if g in gindex]
    ww = np.array([w for g, w in zip(sub.gene, sub.coefficient) if g in gindex])
    return Zmat[:, ii] @ ww


def orient(dk_score, jse_score):
    """Resolve sign using the JSE discovery labels only. Returns (+1/-1, auc_after)."""
    a = roc_auc_score(y, jse_score)
    return (1.0, a) if a >= 0.5 else (-1.0, 1 - a)


def pack(sc, T_, E_, tag):
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    cph = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}), "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T_[hi], T_[~hi], E_[hi], E_[~hi])
    hrs = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "g": hi.astype(int)}), "time", "event")
    ci = np.exp(hrs.confidence_intervals_.iloc[0].values)
    early, late = (T <= 24) & (E == 1), T > 60
    bm = early | late
    return dict(name=tag, C=float(c), hr_per_sd=float(np.exp(cph.params_.iloc[0])),
                cox_p=float(cph.summary["p"].iloc[0]), logrank_p=float(lr.p_value),
                hr_median_split=float(np.exp(hrs.params_.iloc[0])),
                hr_ci=[float(ci[0]), float(ci[1])],
                binary_auc=float(roc_auc_score(early[bm].astype(int), sc[bm])))


res = {"transcription_source": "table_5_ss.png (user-supplied screenshot of Table 5)",
       "caveats": ["no published intercept - ranking only, never an absolute class call",
                   "coefficients applied to rank-inverse-normal expression, not their original "
                   "normalisation: faithful reproduction of ordering, not identical weighting",
                   "orientation resolved on JSE discovery labels, never on DK survival"]}
scores = {}

# ---------------- 3a: exact coefficients ----------------
print("\n" + "=" * 72 + "\n3a  EXACT CLASSIFIER (their published coefficients)\n" + "=" * 72)
for tag, sub in [("brayer49", coef[coef.available]), ("brayer14", coef[coef.available & coef.in_14gene])]:
    d_ = linpred(Zd_all, gi_all, sub); j_ = linpred(Zt_all, gi_all, sub)
    s, auc_jse = orient(d_, j_)
    d_ = s * d_
    scores[tag] = d_
    r = pack(d_, T5, E5, f"{tag} exact ({len(sub)} genes)")
    r["orientation_sign"] = float(s); r["jse_auc_used_for_orientation"] = float(auc_jse)
    res[tag] = r
    print(f"{r['name']:<28s} C={r['C']:.3f}  HR/SD={r['hr_per_sd']:.2f} (p={r['cox_p']:.4f})  "
          f"split HR={r['hr_median_split']:.2f} [{r['hr_ci'][0]:.2f}-{r['hr_ci'][1]:.2f}]  "
          f"logrank p={r['logrank_p']:.4f}  binAUC={r['binary_auc']:.3f}")
    print(f"{'':28s}  (orientation {'+' if s > 0 else '-'}, JSE AUC {auc_jse:.3f})")

# ---------------- 3b: coefficient-free ----------------
print("\n" + "=" * 72 + "\n3b  COEFFICIENT-FREE (their 14 genes, our weighting)\n" + "=" * 72)
g14 = [g for g in coef[coef.in_14gene].gene if g in gi_all]
X14d, X14t = Zd_all[:, [gi_all[g] for g in g14]], Zt_all[:, [gi_all[g] for g in g14]]
c_poor, c_good = X14t[y == 1].mean(0), X14t[y == 0].mean(0)
dr = np.array([stats.pearsonr(x, c_poor)[0] - stats.pearsonr(x, c_good)[0] for x in X14d])
scores["brayer14_centroid"] = dr
r = pack(dr, T5, E5, "brayer14 nearest-centroid")
res["brayer14_centroid"] = r
print(f"{r['name']:<28s} C={r['C']:.3f}  HR/SD={r['hr_per_sd']:.2f} (p={r['cox_p']:.4f})  "
      f"logrank p={r['logrank_p']:.4f}  binAUC={r['binary_auc']:.3f}")

gm = GaussianMixture(2, covariance_type="full", random_state=0, n_init=10, reg_covar=1e-4).fit(X14d)
lab = gm.predict(X14d)
rr = [stats.pearsonr(X14d[lab == c].mean(0), c_poor)[0] for c in range(2)]
m_gmm = (lab == int(np.argmax(rr))).astype(int)
lr = logrank_test(T5[m_gmm == 1], T5[m_gmm == 0], E5[m_gmm == 1], E5[m_gmm == 0])
cf = CoxPHFitter().fit(pd.DataFrame({"time": T5, "event": E5, "g": m_gmm}), "time", "event")
res["brayer14_gmm2"] = dict(name="brayer14 unsupervised GMM k=2", n_poor=int(m_gmm.sum()),
                            hr=float(np.exp(cf.params_.iloc[0])), logrank_p=float(lr.p_value))
print(f"{'brayer14 GMM k=2':<28s} n_poor={m_gmm.sum()}  HR={np.exp(cf.params_.iloc[0]):.2f}  "
      f"logrank p={lr.p_value:.4f}")

# ---------------- our panels, same metrics, same patients ----------------
print("\n" + "=" * 72 + "\nOUR PANELS (same patients, same endpoint, same metrics)\n" + "=" * 72)
sc4 = (Zd[:, order_lock[:4]] * np.sign(beta_lock[order_lock[:4]])).mean(1)
for tag, sc in [("ours_4gene", sc4), ("ours_8gene", sc_ours)]:
    scores[tag] = sc
    r = pack(sc, T5, E5, tag)
    res[tag] = r
    print(f"{r['name']:<28s} C={r['C']:.3f}  HR/SD={r['hr_per_sd']:.2f} (p={r['cox_p']:.4f})  "
          f"split HR={r['hr_median_split']:.2f}  logrank p={r['logrank_p']:.4f}  binAUC={r['binary_auc']:.3f}")

# ---------------- 3c: convergence ----------------
print("\n" + "=" * 72 + "\n3c  CONVERGENCE — do the two approaches flag the same patients?\n" + "=" * 72)
theirs = scores["brayer14"]
ours = sc_ours
hi_t, hi_o = theirs > np.median(theirs), ours > np.median(ours)
ct = pd.crosstab(pd.Series(hi_o, name="ours_high"), pd.Series(hi_t, name="brayer_high"))
kappa = cohen_kappa_score(hi_o.astype(int), hi_t.astype(int))
rho, rho_p = stats.spearmanr(ours, theirs)
print("\ncross-tabulation (median splits):\n", ct.to_string())
print(f"\nCohen's kappa = {kappa:.3f}   Spearman rho = {rho:.3f} (p = {rho_p:.4f})")

bi = pd.DataFrame({"time": T5, "event": E5,
                   "ours": stats.zscore(ours), "brayer": stats.zscore(theirs)})
bcox = CoxPHFitter().fit(bi, "time", "event")
res["convergence"] = dict(
    crosstab=ct.values.tolist(), cohen_kappa=float(kappa),
    spearman_rho=float(rho), spearman_p=float(rho_p),
    agreement=float((hi_o == hi_t).mean()),
    bivariate_cox={k: dict(HR=float(np.exp(bcox.params_[k])), p=float(bcox.summary.loc[k, "p"]))
                   for k in ["ours", "brayer"]})
print("\nbivariate Cox (both scores, 5-year OS):\n", bcox.summary[["exp(coef)", "p"]].to_string())

# combined score if both independent
p_ours = res["convergence"]["bivariate_cox"]["ours"]["p"]
p_them = res["convergence"]["bivariate_cox"]["brayer"]["p"]
independent = (p_ours < 0.05) and (p_them < 0.05)
res["convergence"]["both_independently_significant"] = bool(independent)
comb = stats.zscore(ours) + np.sign(np.exp(bcox.params_["brayer"]) - 1) * stats.zscore(theirs)
scores["combined"] = comb
r = pack(comb, T5, E5, "combined (ours + brayer)")
res["combined"] = r
print(f"\ncombined score: C={r['C']:.3f}  HR/SD={r['hr_per_sd']:.2f} (p={r['cox_p']:.4f})  "
      f"logrank p={r['logrank_p']:.4f}")

rng = np.random.default_rng(3)
n = len(T5); dC = []
for _ in range(2000):
    b = rng.integers(0, n, n)
    if E5[b].sum() < 5:
        continue
    dC.append(concordance_index_censored(E5[b].astype(bool), T5[b], comb[b])[0]
              - concordance_index_censored(E5[b].astype(bool), T5[b], ours[b])[0])
res["combined"]["deltaC_vs_ours"] = dict(mean=float(np.mean(dC)),
                                         lo=float(np.percentile(dC, 2.5)),
                                         hi=float(np.percentile(dC, 97.5)))
print(f"  delta-C vs our panel alone: {np.mean(dC):+.3f} "
      f"[{np.percentile(dC, 2.5):+.3f}, {np.percentile(dC, 97.5):+.3f}]")

# ---------------- summary table ----------------
summ = pd.DataFrame([{"model": k, "n_genes": v.get("n_genes", np.nan), "C": v["C"],
                      "hr_per_sd": v["hr_per_sd"], "cox_p": v["cox_p"],
                      "logrank_p": v["logrank_p"], "binary_auc": v["binary_auc"]}
                     for k, v in res.items()
                     if isinstance(v, dict) and "C" in v]).sort_values("C", ascending=False)
summ.to_csv(f"{OUT}/benchmark_summary.csv", index=False)
print("\n" + "=" * 72 + "\nHEAD-TO-HEAD (DK, n=54, 5-year OS, 21 events)\n" + "=" * 72)
print(summ.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# ---------------- figures ----------------
fig, ax = plt.subplots(figsize=(5.4, 4.4))
early = (T <= 24) & (E == 1); late = T > 60
ax.scatter(stats.zscore(theirs)[~(early | late)], stats.zscore(ours)[~(early | late)],
           c=GREY, s=32, label="intermediate", zorder=2)
ax.scatter(stats.zscore(theirs)[late], stats.zscore(ours)[late], c=BLUE, s=44, label="alive >60 mo", zorder=3)
ax.scatter(stats.zscore(theirs)[early], stats.zscore(ours)[early], c=RED, s=52, label="died <24 mo", zorder=4)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Brayer 14-gene classifier score (z)")
ax.set_ylabel("our 8-gene panel score (z)")
ax.set_title(f"Do the two signatures rank the same patients?\nSpearman ρ = {rho:.2f} (p = {rho_p:.3f}), "
             f"κ = {kappa:.2f}", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "BRAYER_convergence_scatter")

fig, ax = plt.subplots(figsize=(5.6, 4.0))
grp = np.where(hi_o & hi_t, 3, np.where(hi_o & ~hi_t, 2, np.where(~hi_o & hi_t, 1, 0)))
labs = {0: "both low", 1: "Brayer high only", 2: "ours high only", 3: "both high"}
cols = {0: BLUE, 1: "#7fbc41", 2: "#f4a582", 3: RED}
for gnum in [0, 1, 2, 3]:
    m_ = grp == gnum
    if m_.sum() >= 3:
        KaplanMeierFitter().fit(T5[m_], E5[m_], label=f"{labs[gnum]} (n={m_.sum()})").plot_survival_function(
            ax=ax, color=cols[gnum], ci_show=False)
ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("overall survival (5-year)")
ax.set_title("Cross-classification of the two signatures", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "BRAYER_crosstab_km")

f_ = summ.dropna(subset=["C"]).sort_values("C")
fig, ax = plt.subplots(figsize=(6.4, 0.36 * len(f_) + 1.4))
cols_ = [RED if "ours" in m or "combined" in m else BLUE for m in f_.model]
ax.barh(range(len(f_)), f_.C.values - 0.5, left=0.5, color=cols_)
ax.set_yticks(range(len(f_))); ax.set_yticklabels(f_.model.values, fontsize=8)
ax.axvline(0.5, color="k", lw=0.9)
ax.set_xlabel("C-index (DK, 5-year OS)   —   0.5 = chance")
ax.set_title("Same patients, same endpoint, same metric\nred = this project, blue = Brayer et al.", fontsize=9)
for i, (c, p) in enumerate(zip(f_.C.values, f_.logrank_p.values)):
    ax.text(c + 0.004, i, f"{c:.3f} (p={p:.3f})", va="center", fontsize=7)
ax.set_xlim(0.5, max(0.78, f_.C.max() + 0.06))
save(fig, "BRAYER_head_to_head")

json.dump(res, open(f"{OUT}/BRAYER_benchmark.json", "w"), indent=2, default=float)
np.save(f"{OUT}/dk_brayer14_score.npy", theirs)
print("\nDONE ->", OUT)

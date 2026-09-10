"""
Stage 6 - external test of the locked 4-gene panel in the Frerich Oncotarget-2018 ACC cohort.

Status is fixed by work/PRESPEC_stage6_RETROSPECTIVE.md. That file is NOT a pre-specification:
this cohort's outcome label had already been seen when the stage was designed. Read it first.

TEST 1 (primary)   AUC of the 4-gene score for the paper's poor-outcome Group 1 (PO, n = 14)
                   versus Group 2 (n = 52), using the cohort-batch transform.
TEST 2 (secondary) The same contrast under the frozen single-sample transform, and the
                   discovery-locked threshold applied unchanged.

This is a SUBTYPE test. Group 1 was defined by clustering this same expression matrix, and no
survival times exist for these patients. It is the same class of evidence as Stage 5, and it is
NOT independent of Stage 5 - both cohorts draw on the same MD Anderson biorepository (section 4
of the PRESPEC, and TEST 4 below).

Nothing is refitted. Panel, signs, transform, frozen reference and threshold come from
work/release/model_spec.json. Reads only; all outputs go to work/outputs_frerich/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"

ROOT = f"{ACC_DATA_ROOT}"
EXT = f"{ACC_DATA_ROOT}/external_frerich"
OUT = f"{ACC_DATA_ROOT}/model_outputs_frerich"
sys.path.insert(0, HERE)
from score_acc4 import PANEL, score_batch, score_frozen, THRESHOLD  # noqa: E402

RNG = np.random.default_rng(20260827)
B_NULL = 10_000
B_PERM = 10_000
B_BOOT = 5_000

GROUPS = ["PO", "avg", "orange"]
COL = {"PO": "#B03A2E", "avg": "#5D6D7E", "orange": "#E08214"}
CHERRY = "#8A1538"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def boot_auc(y, s, n=B_BOOT):
    out, idx = [], np.arange(len(y))
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        out.append(roc_auc_score(y[b], s[b]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# ---------------------------------------------------------------- load
def load_frerich():
    df = pd.read_csv(f"{EXT}/frerich_hg38_TXcases.csv")
    df.columns = [c.strip() for c in df.columns]
    cols = df.columns[5:].tolist()
    X = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X.index = df["Symbol"].astype(str).values
    empty = [c for c in cols if X[c].sum() == 0]
    X = X.drop(columns=empty)
    X = X.groupby(level=0).sum()
    cpm = np.log2(X.div(X.sum(axis=0), axis=1) * 1e6 + 1).T   # samples x genes
    grp = pd.Series([c.split("_case")[0] for c in cpm.index], index=cpm.index)
    return X.T, cpm, grp, empty


raw, cpm, grp, empty = load_frerich()
y = (grp == "PO").astype(int).values
res = {"cohort": "Frerich 2018 (Oncotarget 9:7341), MD Anderson tissue, UNM sequencing",
       "n_columns_in_file": 68, "dropped_all_zero": empty, "n_analysed": int(len(cpm)),
       "n_PO": int(y.sum()), "n_rest": int((1 - y).sum()),
       "endpoint": "paper Group 1 (poor outcome) vs Group 2 - SUBTYPE label, clustering-derived",
       "prespecified": False, "independent_of_stage5": False}
print(f"n = {len(cpm)}  ({res['n_PO']} PO / {res['n_rest']} rest); dropped {empty}")

# --------------------------------------------------- label verification
myb = cpm["MYB"]
res["label_check_MYB_mean"] = {g: round(float(myb[grp == g].mean()), 2) for g in GROUPS}
res["label_check_MYBL1_mean"] = {g: round(float(cpm["MYBL1"][grp == g].mean()), 2) for g in GROUPS}
res["orange_is_MYB_MYBL1_negative"] = bool(
    (raw.loc[grp == "orange", "MYB"] < 200).all() and (raw.loc[grp == "orange", "MYBL1"] < 200).all())

# ----------------------------------------------------------- TEST 1
s_batch = pd.Series(score_batch(cpm[PANEL]), index=cpm.index)
auc1 = roc_auc_score(y, s_batch.values)
mw = stats.mannwhitneyu(s_batch[y == 1], s_batch[y == 0], alternative="greater")
lo, hi = boot_auc(y, s_batch.values)
perm = np.array([roc_auc_score(RNG.permutation(y), s_batch.values) for _ in range(B_PERM)])
res["test1_batch"] = {
    "AUC": float(auc1), "ci95": [lo, hi],
    "mannwhitney_p_onesided": float(mw.pvalue),
    "label_permutation_p": float((np.sum(perm >= auc1) + 1) / (B_PERM + 1)),
    "mean_score_by_group": {g: round(float(s_batch[grp == g].mean()), 4) for g in GROUPS},
}
print(f"TEST 1  batch AUC {auc1:.3f} [{lo:.3f}-{hi:.3f}]  MW p {mw.pvalue:.3g}  "
      f"perm p {res['test1_batch']['label_permutation_p']:.4f}")

res["per_gene_auc"] = {g: round(float(roc_auc_score(y, cpm[g])), 3) for g in PANEL}
res["gene_detection_zero_samples"] = {g: int((raw[g] == 0).sum()) for g in PANEL}

# secondary: PO vs avg only (orange excluded)
m = (grp != "orange").values
res["test1_PO_vs_avg_only"] = {
    "n": int(m.sum()),
    "AUC": round(float(roc_auc_score(y[m], score_batch(cpm.loc[m, PANEL]))), 3)}

# ----------------------------------------------------------- TEST 2
spec = json.load(open(f"{HERE}/release/model_spec.json"))
ref = pd.DataFrame(spec["frozen_reference"]["values"])[PANEL]
s_froz = pd.Series(score_frozen(cpm[PANEL], ref), index=cpm.index)
auc2 = roc_auc_score(y, s_froz.values)
hi_call = (s_froz.values > THRESHOLD).astype(int)
tp = int(((hi_call == 1) & (y == 1)).sum()); fp = int(((hi_call == 1) & (y == 0)).sum())
fn = int(((hi_call == 0) & (y == 1)).sum()); tn = int(((hi_call == 0) & (y == 0)).sum())
orr, pf = stats.fisher_exact([[tp, fp], [fn, tn]], alternative="greater")
res["test2_frozen"] = {
    "AUC": float(auc2), "threshold": THRESHOLD,
    "tp": tp, "fp": fp, "fn": fn, "tn": tn,
    "sensitivity": round(tp / (tp + fn), 3), "specificity": round(tn / (tn + fp), 3),
    "ppv": round(tp / (tp + fp), 3), "npv": round(tn / (tn + fn), 3),
    "odds_ratio": float(orr), "fisher_p_onesided": float(pf)}
print(f"TEST 2  frozen AUC {auc2:.3f}  sens {res['test2_frozen']['sensitivity']:.2f} "
      f"spec {res['test2_frozen']['specificity']:.2f}  OR {orr:.2f} p {pf:.4g}")

# ------------------------------------------- TEST 3  matched random panels
mu, sd = cpm.mean(), cpm.std()
pool = cpm[cpm.columns[(cpm > 0).mean() >= 0.25]]
mu_p, sd_p = pool.mean(), pool.std()
cands = {}
for g in PANEL:
    d = ((np.log1p(mu_p) - np.log1p(mu[g])) ** 2 / np.var(np.log1p(mu_p))
         + (sd_p - sd[g]) ** 2 / np.var(sd_p))
    cands[g] = d.drop(index=[x for x in PANEL if x in d.index]).nsmallest(400).index.to_numpy()
null = []
while len(null) < B_NULL:
    pick = [RNG.choice(cands[g]) for g in PANEL]
    if len(set(pick)) < 4:
        continue
    null.append(roc_auc_score(y, score_batch(cpm[list(pick)].rename(columns=dict(zip(pick, PANEL))))))
null = np.array(null)
res["test3_matched_nulls"] = {
    "B": B_NULL, "b_at_least_observed": int((null >= auc1).sum()),
    "empirical_p": float((np.sum(null >= auc1) + 1) / (B_NULL + 1)),
    "null_median": float(np.median(null)), "null_p95": float(np.percentile(null, 95)),
    "null_max": float(null.max())}
print(f"TEST 3  matched-null p {res['test3_matched_nulls']['empirical_p']:.4f} "
      f"(median {np.median(null):.3f})")

# ------------------------------------------- TEST 4  overlap with CCR2020
meta = pd.read_csv(f"{EXT}/frerich_case_metadata.csv").set_index("csv_column")
ccl = pd.read_excel(f"{ACC_DATA_ROOT}/external_ccr2020/CCR2020_Clinical.xlsx").set_index("TID")
ccx = pd.read_excel(f"{ACC_DATA_ROOT}/external_ccr2020/ACC_RNAseq.xlsx")
gc = ccx.columns[0]
ccx = (ccx.set_index(ccx[gc].astype(str))
          .drop(columns=[c for c in ccx.columns if ccx[c].dtype == object])
          .groupby(level=0).max())
cclog = np.log2(ccx + 1)
sh = cpm.T.index.intersection(cclog.index)
A, Bm = cpm.T.loc[sh], cclog.loc[sh]
keep = (A.mean(axis=1) > 1) & (Bm.mean(axis=1) > 1)
A, Bm = A[keep], Bm[keep]


def dbl_z(M):
    C = M.sub(M.mean(axis=1), axis=0)
    C = C.sub(C.mean(axis=0), axis=1)
    return C / C.std(axis=0, ddof=0)


Az, Bz = dbl_z(A), dbl_z(Bm)
R = pd.DataFrame((Az.values.T @ Bz.values) / len(Az), index=A.columns, columns=Bm.columns)
pairs = []
for i in R.index:
    a = meta.loc[i]
    if pd.isna(a["age"]):
        continue
    top = R.loc[i].idxmax()
    b = ccl.loc[top]
    if a["sex"][0].upper() == b["Gender"] and abs(a["age"] - b["Age Dx"]) <= 1:
        pairs.append(dict(frerich=i, ccr2020=top, corr=round(float(R.loc[i, top]), 3),
                          sex=b["Gender"], age_frerich=int(a["age"]), age_ccr=int(b["Age Dx"]),
                          acc_subtype=int(b["ACC Subtype"]), ttdeath=str(b["TTDeath"]),
                          status=str(b["Status"])))
obs = len(pairs)
have = [i for i in R.index if pd.notna(meta.loc[i, "age"])]
nulls = []
for _ in range(2000):
    perm_meta = meta.loc[have].sample(frac=1, random_state=int(RNG.integers(1e9)))
    perm_meta.index = have
    c = 0
    for i in have:
        a = perm_meta.loc[i]; b = ccl.loc[R.loc[i].idxmax()]
        if a["sex"][0].upper() == b["Gender"] and abs(a["age"] - b["Age Dx"]) <= 1:
            c += 1
    nulls.append(c)
nulls = np.array(nulls)
res["test4_ccr2020_overlap"] = {
    "concordant_top_matches": obs, "null_mean": float(nulls.mean()),
    "null_p95": float(np.percentile(nulls, 95)), "null_max": int(nulls.max()),
    "p": float((nulls >= obs).mean()), "pairs": pairs}
pd.DataFrame(pairs).to_csv(f"{OUT}/frerich_ccr2020_shared_patients.csv", index=False)
np.savez_compressed(f"{OUT}/frerich_null_draws.npz", matched_panel_auc=null, overlap_counts=nulls)
print(f"TEST 4  {obs} metadata-concordant top matches vs null mean {nulls.mean():.2f} "
      f"(p < {max((nulls >= obs).mean(), 1/2000):.4f})")

# label concordance on the shared patients
pp = pd.DataFrame(pairs)
pp["po"] = pp.frerich.str.startswith("PO").astype(int)
a_ = int(((pp.po == 1) & (pp.acc_subtype == 1)).sum()); b_ = int(((pp.po == 1) & (pp.acc_subtype == 2)).sum())
c_ = int(((pp.po == 0) & (pp.acc_subtype == 1)).sum()); d_ = int(((pp.po == 0) & (pp.acc_subtype == 2)).sum())
_, pfc = stats.fisher_exact([[a_, b_], [c_, d_]])
res["test4_label_concordance"] = {"PO_ACCI": a_, "PO_ACCII": b_, "rest_ACCI": c_,
                                  "rest_ACCII": d_, "fisher_p": float(pfc)}

# ------------------------------------------- confounding descriptives
res["score_correlations"] = {}
for g in ["MYB", "MYBL1", "SOX4", "MKI67", "TOP2A"]:
    r = stats.spearmanr(s_batch.values, cpm[g].values)
    res["score_correlations"][g] = {"rho": round(float(r.statistic), 3),
                                    "p": float(r.pvalue),
                                    "own_auc": round(float(roc_auc_score(y, cpm[g])), 3)}

# ------------------------------------------------------------ save scores
pd.DataFrame({"sample": cpm.index, "group": grp.values, "y_poor_outcome": y,
              "score_batch": s_batch.values, "score_frozen": s_froz.values,
              "high_risk_locked_threshold": hi_call}).to_csv(
    f"{OUT}/frerich_scores.csv", index=False)
json.dump(res, open(f"{OUT}/FRERICH_results.json", "w"), indent=2)

# ================================================================== FIGURE
fig = plt.figure(figsize=(13.5, 8.4))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.30)

# (A) score by group
ax = fig.add_subplot(gs[0, 0])
for k, g in enumerate(GROUPS):
    v = s_batch[grp == g].values
    ax.scatter(np.full(len(v), k) + RNG.normal(0, 0.07, len(v)), v, s=26,
               color=COL[g], alpha=.85, edgecolor="white", linewidth=.5, zorder=3)
    ax.hlines(np.median(v), k - .28, k + .28, color="black", lw=2, zorder=4)
ax.set_xticks(range(3))
ax.set_xticklabels([f"Group 1\npoor outcome\n(n={(grp=='PO').sum()})",
                    f"Group 2\n$MYB$/$MYBL1$+\n(n={(grp=='avg').sum()})",
                    f"Group 2\nneither\n(n={(grp=='orange').sum()})"], fontsize=8.5)
ax.set_ylabel("locked four-gene score")
ax.set_title("A   Score by published group", loc="left", fontweight="bold", fontsize=10)
ax.axhline(0, color="#cccccc", lw=.8, zorder=1)
ax.text(.98, .97, f"Mann–Whitney\n$P$ = {mw.pvalue:.2g}", transform=ax.transAxes,
        ha="right", va="top", fontsize=8.5, style="italic")

# (B) ROC
ax = fig.add_subplot(gs[0, 1])
fpr, tpr, _ = roc_curve(y, s_batch.values)
ax.plot(fpr, tpr, color=CHERRY, lw=2.4, label=f"batch score  AUC {auc1:.3f}")
fpr2, tpr2, _ = roc_curve(y, s_froz.values)
ax.plot(fpr2, tpr2, color="#2E86AB", lw=1.8, ls="--", label=f"frozen score  AUC {auc2:.3f}")
ax.plot([0, 1], [0, 1], color="#bbbbbb", lw=1, ls=":")
ax.set_xlabel("1 − specificity"); ax.set_ylabel("sensitivity")
ax.set_title("B   Discrimination of Group 1", loc="left", fontweight="bold", fontsize=10)
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.text(.03, .95, f"95% CI {lo:.3f}–{hi:.3f}", transform=ax.transAxes, fontsize=8, va="top")

# (C) matched nulls
ax = fig.add_subplot(gs[0, 2])
ax.hist(null, bins=44, color="#B8C4CE", edgecolor="white", linewidth=.4)
ax.axvline(auc1, color=CHERRY, lw=2.4)
ax.text(auc1, ax.get_ylim()[1] * .93, f"  locked panel\n  AUC {auc1:.3f}",
        color=CHERRY, fontsize=8.5, fontweight="bold", va="top")
ax.set_xlabel("AUC of matched random four-gene panel"); ax.set_ylabel("count")
ax.set_title("C   Expression- and variance-matched null", loc="left", fontweight="bold", fontsize=10)
ax.text(.03, .72, f"$B$ = {B_NULL:,}\n$b$ = {res['test3_matched_nulls']['b_at_least_observed']}\n"
                  f"$P$ = {res['test3_matched_nulls']['empirical_p']:.4f}",
        transform=ax.transAxes, fontsize=8.5, va="top")

# (D) locked threshold
ax = fig.add_subplot(gs[1, 0])
bars = [[tn, fp], [fn, tp]]
ax.bar([0, 1], [tn, fn], color="#5D6D7E", edgecolor="white", label="called low risk")
ax.bar([0, 1], [fp, tp], bottom=[tn, fn], color=CHERRY, edgecolor="white", label="called high risk")
for k, (a_, b_) in enumerate(zip([tn, fn], [fp, tp])):
    ax.text(k, a_ / 2, str(a_), ha="center", va="center", color="white", fontweight="bold")
    ax.text(k, a_ + b_ / 2, str(b_), ha="center", va="center", color="white", fontweight="bold")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Group 2\n(n=52)", "Group 1\n(n=14)"], fontsize=9)
ax.set_ylabel("patients")
ax.set_title("D   Discovery-locked threshold", loc="left", fontweight="bold", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.text(.97, .62, f"sens {res['test2_frozen']['sensitivity']:.2f}\n"
                  f"spec {res['test2_frozen']['specificity']:.2f}\n"
                  f"OR {orr:.1f}\n$P$ = {pf:.4f}", transform=ax.transAxes,
        ha="right", fontsize=8.5, va="top")

# (E) per-gene
ax = fig.add_subplot(gs[1, 1])
vals = [res["per_gene_auc"][g] for g in PANEL]
ax.barh(range(4), vals, color="#8FA6B8", edgecolor="white", height=.62)
ax.axvline(auc1, color=CHERRY, lw=2, label=f"panel {auc1:.3f}")
ax.axvline(.5, color="#bbbbbb", lw=1, ls=":")
ax.set_yticks(range(4)); ax.set_yticklabels(PANEL, style="italic")
ax.set_xlim(.4, .9); ax.set_xlabel("AUC, Group 1 vs Group 2")
ax.set_title("E   Single-gene discrimination", loc="left", fontweight="bold", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper right")
for k, v in enumerate(vals):
    ax.text(v + .006, k, f"{v:.3f}", va="center", fontsize=8)

# (F) overlap with CCR2020
ax = fig.add_subplot(gs[1, 2])
ax.hist(nulls, bins=range(0, max(nulls.max(), obs) + 2), color="#B8C4CE",
        edgecolor="white", linewidth=.4, align="left")
ax.axvline(obs, color=CHERRY, lw=2.4)
ax.text(obs, ax.get_ylim()[1] * .93, f" observed\n {obs} shared",
        color=CHERRY, fontsize=8.5, fontweight="bold", va="top")
ax.set_xlabel("metadata-concordant top expression matches")
ax.set_ylabel("count")
ax.set_title("F   Patient overlap with CCR2020", loc="left", fontweight="bold", fontsize=10)
ax.text(.55, .60, f"null mean {nulls.mean():.1f}\n95th pct {np.percentile(nulls,95):.0f}\n"
                  f"$P$ < {max((nulls>=obs).mean(),1/2000):.4f}",
        transform=ax.transAxes, fontsize=8.5, va="top")

for a in fig.get_axes():
    a.spines[["top", "right"]].set_visible(False)
save(fig, "FIG_FRERICH_external")
print("wrote", OUT)

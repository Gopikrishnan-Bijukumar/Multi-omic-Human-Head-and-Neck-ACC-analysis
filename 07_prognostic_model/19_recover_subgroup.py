"""STAGE 3d: can OUR candidate biomarker panel identify THEIR aggressive subgroup?

Different question from 18_/18b_. Those asked whether the two scores agree. This one treats the
Brayer subgroup as an answer key and asks how well our 4-gene panel recovers it.

  Q1  Recovery      - AUC of our panel score for predicting membership of their subgroup,
                      across every subgroup size, with a null from 2000 random 4-gene panels.
  Q2  Symmetry      - the reverse: how well does their classifier recover OUR high-risk group?
  Q3  Discordance   - when the two disagree, who actually died? This is the one that matters:
                      it says which signature is right where they differ.

Q1 and Q2 use EXPRESSION ONLY - no survival involved, so no multiplicity concern there.
Q3 uses survival and is descriptive/exploratory; labelled as such, and it rests on few patients.
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
from sklearn.metrics import roc_auc_score, roc_curve
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_recovery"
import os; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY, GREEN, ORANGE = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837", "#e08214"


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw, clin = D["dk_raw"], D["tr_raw"], D["clin"]
keep = (dk_raw.mean(0).values > 1.0) & ((dk_raw > 0).mean(0).values > 0.6) & (tr_raw.mean(0).values > 1.0)
genes = [g for g, k in zip(genes_all, keep) if k]
Zd = d3.rint(dk_raw[genes]).values
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
n = len(T)

LOCK = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta_lock = pd.Series(np.array(LOCK["beta"]), index=LOCK["genes"]).reindex(genes).values
order_lock = np.argsort(-np.abs(beta_lock))
idx4, idx8 = order_lock[:4], order_lock[:8]
panel4 = [genes[i] for i in idx4]


def pscore(idx):
    return (Zd[:, idx] * np.sign(beta_lock[idx])).mean(1)


ours4, ours8 = pscore(idx4), pscore(idx8)
theirs = np.load(f"{ACC_DATA_ROOT}/model_outputs_brayer/dk_brayer14_score.npy")
print(f"DK n={n}, events(5y)={E5.sum()}")
print(f"our candidate panel: {panel4}")

res = {"panel4": panel4, "n_dk": int(n)}

# =====================================================================
# Q1 — how well does our panel recover their subgroup, at every size?
# =====================================================================
SIZES = [4, 6, 8, 10, 12, 14, 16, 20, 24, 27]
rng = np.random.default_rng(23)

rows = []
for k in SIZES:
    thr = np.sort(theirs)[-k]
    their_grp = (theirs >= thr).astype(int)
    if their_grp.sum() != k:                       # ties
        their_grp = np.zeros(n, int); their_grp[np.argsort(-theirs)[:k]] = 1

    auc4 = roc_auc_score(their_grp, ours4)
    auc8 = roc_auc_score(their_grp, ours8)

    # matched-size overlap: our own top-k vs their top-k
    our_top = np.zeros(n, bool); our_top[np.argsort(-ours4)[:k]] = True
    ov = int((our_top & their_grp.astype(bool)).sum())
    hg = stats.hypergeom.sf(ov - 1, n, k, k)
    sens = ov / k
    jac = ov / len(set(np.where(our_top)[0]) | set(np.where(their_grp)[0]))

    # null: 2000 random 4-gene sign-weighted panels
    null = np.array([roc_auc_score(
        their_grp, (Zd[:, rng.choice(len(genes), 4, replace=False)] *
                    rng.choice([-1.0, 1.0], 4)).mean(1)) for _ in range(2000)])
    p_null = (np.sum(null >= auc4) + 1) / 2001

    rows.append(dict(subgroup_size=k, frac=k / n, auc_panel4=auc4, auc_panel8=auc8,
                     overlap=ov, sensitivity=sens, jaccard=jac,
                     hypergeom_p=float(hg), null_p=float(p_null),
                     null_mean=float(null.mean()), null_p95=float(np.percentile(null, 95))))
rec = pd.DataFrame(rows)
rec.to_csv(f"{OUT}/recovery_by_subgroup_size.csv", index=False)
print("\nQ1 — recovering the Brayer subgroup with our 4-gene panel (expression only):")
print(rec[["subgroup_size", "auc_panel4", "auc_panel8", "overlap", "sensitivity",
           "hypergeom_p", "null_p"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
res["Q1_recovery"] = rec.to_dict("records")

# =====================================================================
# Q2 — symmetry: does their classifier recover OUR high-risk group?
# =====================================================================
sym = []
for k in SIZES:
    our_grp = np.zeros(n, int); our_grp[np.argsort(-ours4)[:k]] = 1
    sym.append(dict(subgroup_size=k, auc_theirs_recovering_ours=roc_auc_score(our_grp, theirs)))
sym = pd.DataFrame(sym)
sym.to_csv(f"{OUT}/symmetry.csv", index=False)
res["Q2_symmetry"] = sym.to_dict("records")
print("\nQ2 — reverse direction (their classifier recovering our high-risk group):")
print(sym.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

# =====================================================================
# Q3 — DISCORDANCE: when they disagree, who dies?
# =====================================================================
print("\nQ3 — discordance analysis (descriptive; small numbers)")
disc_all = {}
for k in [6, 11, 14, 27]:
    o = np.zeros(n, bool); o[np.argsort(-ours4)[:k]] = True
    t = np.zeros(n, bool); t[np.argsort(-theirs)[:k]] = True
    cells = {"both": o & t, "ours_only": o & ~t, "theirs_only": ~o & t, "neither": ~o & ~t}
    d = {}
    for name, m in cells.items():
        if m.sum() == 0:
            continue
        km = KaplanMeierFitter().fit(T[m], E[m])
        d[name] = dict(n=int(m.sum()),
                       deaths_5yr=int(E5[m].sum()),
                       death_rate_5yr=float(E5[m].mean()),
                       median_os=float(km.median_survival_time_))
    disc_all[f"top{k}"] = d
    print(f"\n  top-{k} groups:")
    for name, v in d.items():
        print(f"    {name:<12s} n={v['n']:>2d}  deaths<60mo={v['deaths_5yr']:>2d} "
              f"({v['death_rate_5yr']*100:>3.0f}%)  median OS={v['median_os']:.0f} mo")
res["Q3_discordance"] = disc_all

# formal test on the discordant cells at the most informative size
k = 11
o = np.zeros(n, bool); o[np.argsort(-ours4)[:k]] = True
t = np.zeros(n, bool); t[np.argsort(-theirs)[:k]] = True
oo, to = o & ~t, ~o & t
if oo.sum() >= 3 and to.sum() >= 3:
    lr = logrank_test(T5[oo], T5[to], E5[oo], E5[to])
    res["Q3_discordant_headtohead"] = dict(
        k=k, n_ours_only=int(oo.sum()), n_theirs_only=int(to.sum()),
        deaths_ours_only=int(E5[oo].sum()), deaths_theirs_only=int(E5[to].sum()),
        logrank_p=float(lr.p_value))
    print(f"\n  head-to-head on discordant patients (top-{k}): "
          f"ours-only {E5[oo].sum()}/{oo.sum()} died vs theirs-only {E5[to].sum()}/{to.sum()} died, "
          f"log-rank p={lr.p_value:.4f}")

# =====================================================================
# figures
# =====================================================================
fig, ax = plt.subplots(figsize=(5.6, 3.9))
ax.fill_between(rec.subgroup_size, rec.null_mean, rec.null_p95, color=GREY, alpha=.3,
                label="random 4-gene panels (mean–95th pct)")
ax.plot(rec.subgroup_size, rec.auc_panel4, "o-", color=RED, ms=5, label="our 4-gene panel")
ax.plot(rec.subgroup_size, rec.auc_panel8, "s--", color=ORANGE, ms=4, label="our 8-gene panel")
ax.axhline(0.5, color="k", lw=.8, ls=":")
ax.set_xlabel("size of the Brayer aggressive subgroup (patients)")
ax.set_ylabel("AUC recovering that subgroup")
ax.set_ylim(0.35, 1.0)
ax.set_title("Can our panel identify their aggressive subgroup?\n(expression only — no survival used)", fontsize=9)
ax.legend(fontsize=8, frameon=False, loc="lower right")
save(fig, "RECOVERY_auc_by_size")

# ROC at their best-performing subgroup size (n=6)
k = 6
their_grp = np.zeros(n, int); their_grp[np.argsort(-theirs)[:k]] = 1
fig, ax = plt.subplots(figsize=(4.3, 4.1))
for sc, lab, col in [(ours4, "our 4-gene panel", RED), (ours8, "our 8-gene panel", ORANGE)]:
    fpr, tpr, _ = roc_curve(their_grp, sc)
    ax.plot(fpr, tpr, color=col, lw=1.8, label=f"{lab} (AUC {roc_auc_score(their_grp, sc):.3f})")
ax.plot([0, 1], [0, 1], "k:", lw=.9)
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
ax.set_title(f"Recovering the Brayer top-{k} subgroup", fontsize=9)
ax.legend(fontsize=8, frameon=False, loc="lower right")
save(fig, "RECOVERY_roc")

# discordance KM at k=11
k = 11
o = np.zeros(n, bool); o[np.argsort(-ours4)[:k]] = True
t = np.zeros(n, bool); t[np.argsort(-theirs)[:k]] = True
fig, ax = plt.subplots(figsize=(5.6, 4.0))
for name, m, col in [("both flag (n=%d)" % (o & t).sum(), o & t, RED),
                     ("ours only (n=%d)" % (o & ~t).sum(), o & ~t, ORANGE),
                     ("theirs only (n=%d)" % (~o & t).sum(), ~o & t, GREEN),
                     ("neither (n=%d)" % (~o & ~t).sum(), ~o & ~t, BLUE)]:
    if m.sum() >= 2:
        KaplanMeierFitter().fit(T5[m], E5[m], label=name).plot_survival_function(
            ax=ax, color=col, ci_show=False)
ax.set_ylim(0, 1.02); ax.set_xlabel("months"); ax.set_ylabel("overall survival (5-year)")
ax.set_title(f"When the two signatures disagree, who dies?\n(top-{k} by each; descriptive)", fontsize=9)
ax.legend(fontsize=8, frameon=False)
save(fig, "RECOVERY_discordance_km")

json.dump(res, open(f"{OUT}/RECOVERY_results.json", "w"), indent=2, default=float)
print("\nDONE ->", OUT)

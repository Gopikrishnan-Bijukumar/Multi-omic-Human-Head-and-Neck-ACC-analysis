"""Stage 36c (v2 of stage 31c) - manuscript figures rebuilt from outputs_metrics_v2.

Identical to 31c except for the output directory.  See 36_unified_metrics_v2.py for why
v2 exists.  work/outputs_metrics/ figures are left untouched.

Original stage 31c docstring follows.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)


"""STAGE 7c - paper figures for the unified metrics.

Split from 31_unified_metrics.py so figures can be iterated without re-running the 2,000-resample
bootstraps. Reads work/outputs_metrics/UNIFIED_metrics.json and the score vectors.

Figures
  M1  both metrics, same patients   (A) paired C / tdAUC forest   (B) IPCW ROC at 60 months
  M2  endpoint-blocked cross-cohort summary - the figure that stops 0.884 and 0.671 being read
      as comparable
  M3  discrimination over follow-up - tdAUC(t) with Harrell C as reference
  S1  metric agreement - C vs tdAUC, Harrell vs Uno

Every figure is written twice: .png at 600 dpi and .svg with svg.fonttype "none". All output goes
through save() - no bare fig.savefig anywhere.

Palette: Variant B (the 23-30 review-response vintage), matching the neighbouring DELONG_*,
CCR2020_* and THRESHOLD_* figures.

Terminology: the discovery cohort is the "internal cohort" in every label. The on-disk column
`score_frozen_jse` is mapped to "frozen internal reference" at plot time; the file is untouched.

Reads only; writes to work/outputs_metrics/.
"""
import importlib.util as _il
import json
import os

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from sklearn.metrics import roc_auc_score
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_metrics_v2"   # v2: reads/writes the regenerated metrics only

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

# Variant B palette, extended to six series
COL = {
    "ours_4gene":        "#B4436C",
    "ours_4gene_frozen": "#8F3556",
    "ours_8gene":        "#C9788F",
    "composition":       "#7D8F3C",
    "brayer14":          "#3C6E8F",
    "brayer49":          "#6E8FA8",
}
GREY = "#9a9a9a"
# CCR2020 uses its own score keys; map them onto the same colours
COL.update({"score4": COL["ours_4gene"], "score8": COL["ours_8gene"],
            "myc_tp63": GREY, "brayer14_ccr": COL["brayer14"]})


def save(fig, name):
    """Every figure in both formats; PNG at 600 dpi, SVG with live text."""
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def despine(ax):
    ax.spines[["top", "right"]].set_visible(False)


def panel_label(ax, letter):
    """A/B/C lettering - a new convention for this project, introduced for the paper figures."""
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="left")


def _load(name, path):
    s = _il.spec_from_file_location(name, path)
    m = _il.module_from_spec(s)
    s.loader.exec_module(m)
    return m


d3 = _load("d3", f"{HERE}/03_data.py")

M = json.load(open(f"{OUT}/UNIFIED_metrics.json"))
dk = pd.DataFrame(M["dk"])
ccr = pd.DataFrame(M["ccr2020"])
bindf = pd.DataFrame(M["binary"])
TIMES = np.array(M["tdauc_times_months"], float)
CURVES = {k: np.array(v, float) for k, v in M["tdauc_curves"].items()}

# ---- rebuild the DK score vectors (cheap; no bootstraps)
D = d3.load("all")
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)

s4 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
s8 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel8_score.npy")
sb14 = np.load(f"{ACC_DATA_ROOT}/model_outputs_brayer/dk_brayer14_score.npy")
thr = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_threshold/dk_scores_by_procedure.csv")
s_frozen = thr["score_frozen_jse"].values

coef = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_brayer/brayer_table5_coefs.csv")
Zd_all, Zt_all = d3.rint(dk_raw).values, d3.rint(tr_raw).values
gi_all = {g: i for i, g in enumerate(genes_all)}
sub = coef[coef.gene.map(lambda g: g in gi_all)]
ii = [gi_all[g] for g in sub.gene]
ww = sub.coefficient.values
sgn = 1.0 if roc_auc_score(y, Zt_all[:, ii] @ ww) >= 0.5 else -1.0
sb49 = sgn * (Zd_all[:, ii] @ ww)

dk_comp = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_deconv/dk_estimated_composition.csv", index_col=0)
TRIAD = {"Epithelial_Tumor": +1.0, "Dividing": +1.0, "Myoepithelial_Tumor": -1.0}
zc = (dk_comp - dk_comp.mean(0)) / (dk_comp.std(0) + 1e-12)
s_comp = sum(w * zc[c] for c, w in TRIAD.items()).values

SC = {"ours_4gene": s4, "ours_4gene_frozen": s_frozen, "ours_8gene": s8,
      "composition": s_comp, "brayer14": sb14, "brayer49": sb49}
ORDER = ["ours_4gene", "ours_4gene_frozen", "ours_8gene", "composition", "brayer14", "brayer49"]
LAB = dict(zip(dk.score, dk.label))


# ---------------------------------------------------------------- IPCW ROC at a horizon
def ipcw_roc(T_, E_, sc, t):
    """Cumulative/dynamic ROC at horizon t, IPCW-weighted (Uno et al.).

    cases    = died by t                 weight 1 / G(T_i-)
    controls = still at risk after t     weight 1 / G(t)
    Patients censored before t contribute to neither, which is exactly what the weights correct
    for. Returns (fpr, tpr, auc_trapezoid).
    """
    km = KaplanMeierFitter().fit(T_, 1 - E_)          # censoring distribution
    def G(x):
        return np.clip(km.predict(x), 1e-8, None)

    case = (T_ <= t) & (E_ == 1)
    ctrl = T_ > t
    wc = 1.0 / np.asarray(G(T_[case]), float)
    wk = np.full(ctrl.sum(), 1.0 / float(G(t)))

    thrs = np.unique(np.concatenate([sc, [sc.min() - 1, sc.max() + 1]]))[::-1]
    tpr = np.array([(wc * (sc[case] > th)).sum() / wc.sum() for th in thrs])
    fpr = np.array([(wk * (sc[ctrl] > th)).sum() / wk.sum() for th in thrs])
    o = np.argsort(fpr)
    return fpr[o], tpr[o], float(np.trapezoid(tpr[o], fpr[o]))


# ================================================================ M1
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6))

# (A) paired forest: C(5yr) and tdAUC@60 on one axis
a = ax[0]
yy = np.arange(len(ORDER))[::-1]
for i, k in enumerate(ORDER):
    r = dk[dk.score == k].iloc[0]
    c = COL[k]
    a.plot([r.C_5yr_ci_lo, r.C_5yr_ci_hi], [yy[i] + .13] * 2, color=c, lw=1.5, alpha=.85, zorder=2)
    a.plot(r.C_5yr, yy[i] + .13, "o", color=c, ms=7, zorder=3)
    a.plot([r.tdAUC_60mo_ci_lo, r.tdAUC_60mo_ci_hi], [yy[i] - .13] * 2, color=c, lw=1.5,
           alpha=.55, zorder=2)
    a.plot(r.tdAUC_60mo, yy[i] - .13, "s", color="white", mec=c, mew=1.6, ms=7, zorder=3)
a.axvline(.5, color=GREY, ls=":", lw=1)
a.set_yticks(yy)
a.set_yticklabels([LAB[k] for k in ORDER], fontsize=8.5)
a.set_xlim(.32, 1.0)
a.set_xlabel("discrimination (0.5 = chance)")
a.set_title("Same 54 patients, same endpoint, two metrics", fontsize=9)
a.legend(handles=[Line2D([], [], marker="o", color=GREY, ls="none", ms=7,
                         label="Harrell C, 5-year OS"),
                  Line2D([], [], marker="s", color=GREY, mfc="white", ls="none", ms=7,
                         label="time-dependent AUC @ 60 mo")],
         fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.13), ncol=2)
despine(a)
panel_label(a, "A")

# (B) IPCW ROC at 60 months - all 54 patients, unlike the n=44 binary ROC in DELONG_*
b = ax[1]
for k in ORDER:
    fpr, tpr, auc_t = ipcw_roc(T, E, SC[k], 60.0)
    r = dk[dk.score == k].iloc[0]
    b.step(fpr, tpr, where="post", color=COL[k], lw=1.7,
           label=f"{LAB[k]}  {r.tdAUC_60mo:.3f}")
b.plot([0, 1], [0, 1], ls=":", color=GREY, lw=1)
b.set_xlabel("1 - specificity")
b.set_ylabel("sensitivity")
b.set_xlim(-.01, 1.01)
b.set_ylim(-.01, 1.01)
b.set_title("Time-dependent ROC at 60 months (IPCW, n = 54)", fontsize=9)
b.legend(fontsize=7.5, frameon=False, loc="lower right")
despine(b)
panel_label(b, "B")

fig.suptitle("DK external validation: the C-index and the AUC are measuring the same thing",
             fontsize=10, y=1.02)
save(fig, "M1_dk_both_metrics")
print("M1 written")

# ================================================================ M2
BLOCKS = [
    ("DK survival  (5-year OS, n = 54)", "survival"),
    ("DK early death vs long survival  (n = 44)", "binary_dk"),
    ("DK Brayer published poor-survival group  (n = 62)", "binary_bg"),
    ("CCR2020 ACC-I vs ACC-II subtype  (n = 54)", "binary_ccr"),
    ("CCR2020 time to death, decedents only  (n = 34)\n"
     "     conditioned on death - NOT a population C, not comparable to DK",
     "ccr_dec"),
    ("Internal cohort, nested CV  (n = 20)", "internal"),
]

rows = []
for k in ORDER:
    r = dk[dk.score == k].iloc[0]
    rows.append(("survival", LAB[k], r.C_5yr, r.C_5yr_ci_lo, r.C_5yr_ci_hi, "C", COL[k], False))
    rows.append(("survival", LAB[k], r.tdAUC_60mo, r.tdAUC_60mo_ci_lo, r.tdAUC_60mo_ci_hi,
                 "AUC", COL[k], False))
for _, r in bindf[bindf.block == "DK early vs late"].iterrows():
    rows.append(("binary_dk", r.label, r.auc, r.ci_lo, r.ci_hi, "AUC",
                 COL.get(r.score, GREY), bool(r.in_sample)))
for _, r in bindf[bindf.block == "DK Brayer published group"].iterrows():
    rows.append(("binary_bg", r.label, r.auc, r.ci_lo, r.ci_hi, "AUC", COL["ours_4gene"], False))
for _, r in bindf[bindf.block == "CCR2020 subtype"].iterrows():
    rows.append(("binary_ccr", r.label, r.auc, r.ci_lo, r.ci_hi, "AUC",
                 COL.get(r.score, GREY) if not r.in_sample else GREY, bool(r.in_sample)))
for _, r in ccr[ccr.score == "score4"].iterrows():
    rows.append(("ccr_dec", "4-gene panel", r.decedent_C, r.decedent_C_ci_lo, r.decedent_C_ci_hi,
                 "C", COL["ours_4gene"], False))
    rows.append(("ccr_dec", "4-gene panel", r.decedent_auc_lt5y, r.decedent_auc_ci_lo,
                 r.decedent_auc_ci_hi, "AUC", COL["ours_4gene"], False))
for _, r in bindf[bindf.block == "Internal cohort (nested CV)"].iterrows():
    rows.append(("internal", r.label, r.auc, r.ci_lo, r.ci_hi, "AUC", COL["ours_4gene"], False))

R = pd.DataFrame(rows, columns=["blk", "label", "v", "lo", "hi", "metric", "col", "ins"])

present = [(t, b) for t, b in BLOCKS if (R.blk == b).any()]
h = 0.30 * len(R) + 0.62 * len(present) + 1.5
fig, a = plt.subplots(figsize=(8.6, h))
ypos, ylabels, seps = [], [], []
cur = 0.0
for title, b in present:
    seps.append((cur + 0.55, title))
    cur -= 0.55
    for _, r in R[R.blk == b].iterrows():
        ypos.append(cur)
        ylabels.append(f"{r.label}   [{r.metric}]" + ("  in-sample" if r.ins else ""))
        cur -= 1.0
    cur -= 0.30
for i, (_, r) in enumerate(R.iterrows()):
    yv = ypos[i]
    fill = "white" if (r.metric == "AUC" or r.ins) else r.col
    a.plot([r.lo, r.hi], [yv, yv], color=r.col, lw=1.6, alpha=.45 if r.ins else .85, zorder=2)
    a.plot(r.v, yv, "s" if r.metric == "AUC" else "o", color=fill, mec=r.col, mew=1.6, ms=7.5,
           alpha=.55 if r.ins else 1.0, zorder=3)
    a.text(min(r.hi, 1.02) + .012, yv, f"{r.v:.3f}", va="center", fontsize=7.4, color=r.col,
           alpha=.55 if r.ins else 1.0)
a.axvline(.5, color=GREY, ls=":", lw=1, zorder=1)
for yv, title in seps:
    a.text(.305, yv, title, fontsize=8.4, fontweight="bold", va="center", ha="left")
a.set_yticks(ypos)
a.set_yticklabels(ylabels, fontsize=7.4)
a.set_ylim(cur + .6, seps[0][0] + .7)
a.set_xlim(.30, 1.12)
a.set_xticks(np.arange(.3, 1.05, .1))
a.set_xlabel("discrimination (0.5 = chance)")
a.set_title("Every validation, blocked by endpoint\n"
            "values are comparable WITHIN a block, never across blocks", fontsize=9.5)
a.legend(handles=[Line2D([], [], marker="o", color=GREY, ls="none", ms=7, label="C-index"),
                  Line2D([], [], marker="s", color=GREY, mfc="white", ls="none", ms=7,
                         label="AUC"),
                  Line2D([], [], marker="s", color=GREY, mfc="white", ls="none", ms=7,
                         alpha=.55, label="in-sample comparator")],
         fontsize=7.8, frameon=False, loc="upper right", bbox_to_anchor=(1.0, .995))
despine(a)
save(fig, "M2_endpoint_blocked_summary")
print("M2 written")

# ================================================================ M3
fig, a = plt.subplots(figsize=(6.6, 4.4))
for k in ORDER:
    r = dk[dk.score == k].iloc[0]
    a.plot(TIMES, CURVES[k], color=COL[k], lw=1.8, label=f"{LAB[k]}")
    a.axhline(r.C_5yr, color=COL[k], ls="--", lw=.8, alpha=.45)
a.axhline(.5, color=GREY, ls=":", lw=1)
a.set_xlabel("months since surgery")
a.set_ylabel("time-dependent AUC")
a.set_ylim(.35, .95)
a.set_xlim(TIMES[0], TIMES[-1])
a.set_title("Discrimination across follow-up\n"
            "solid = tdAUC(t);  dashed = that score's Harrell C over 5 years", fontsize=9)
a.legend(fontsize=7.6, frameon=False, loc="upper right", ncol=2)
despine(a)
save(fig, "M3_tdauc_over_followup")
print("M3 written")

# ================================================================ S1
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.3))

# (A) do the two metrics rank the six scores the same way?
a = ax[0]
for k in ORDER:
    r = dk[dk.score == k].iloc[0]
    a.plot(r.C_5yr, r.tdAUC_60mo, "o", color=COL[k], ms=10, label=LAB[k])
lim = [.46, .76]
a.plot(lim, lim, ls=":", color=GREY, lw=1)
a.set_xlim(*lim)
a.set_ylim(*lim)
a.set_xlabel("Harrell C, 5-year OS")
a.set_ylabel("time-dependent AUC @ 60 mo")
rho = np.corrcoef(dk.C_5yr, dk.tdAUC_60mo)[0, 1]
a.set_title(f"The two metrics agree  (r = {rho:.3f})", fontsize=9)
a.legend(fontsize=7, frameon=False, loc="upper left")
despine(a)
panel_label(a, "A")

# (B) why Harrell and Uno coincide exactly here.
# NOT a null result: DK has ZERO censoring events before 60 months - every patient either died
# within 5 years (21) or was followed beyond 5 years (33). The censoring distribution G(t) is
# therefore 1.0 throughout the window, every IPCW weight is 1, and Uno's C reduces to Harrell's
# C algebraically. The useful consequence is that the 5-year C-index carries no censoring bias.
b = ax[1]
km = KaplanMeierFitter().fit(T, 1 - E)
tt = np.linspace(0, float(T.max()), 400)
b.step(tt, np.asarray(km.predict(tt), float), where="post", color="#3C6E8F", lw=2)
b.axvline(60, color="#B4436C", ls="--", lw=1.4)
b.text(61, .18, "5-year horizon", color="#B4436C", fontsize=8, rotation=90, va="bottom")
b.axvspan(0, 60, color="#3C6E8F", alpha=.08)
b.text(30, .45, "G(t) = 1.000\nthroughout", ha="center", fontsize=8.5, color="#3C6E8F")
n_cens_early = int(((T <= 60) & (E == 0)).sum())
b.set_xlabel("months since surgery")
b.set_ylabel("G(t):  probability of remaining uncensored")
b.set_ylim(0, 1.04)
b.set_xlim(0, float(T.max()))
b.set_title(f"No censoring inside the 5-year window\n"
            f"{n_cens_early} of {len(T)} censored before 60 mo -> every IPCW weight = 1,\n"
            f"so Uno C = Harrell C exactly (max diff "
            f"{float(np.abs(dk.C_5yr - dk.C_uno_tau60).max()):.1e})", fontsize=8.4)
despine(b)
panel_label(b, "B")
fig.suptitle("The 5-year C-index in DK is free of censoring bias", fontsize=10, y=1.03)
save(fig, "S1_metric_agreement")
print("S1 written")

print(f"\nfigures written -> {OUT}")

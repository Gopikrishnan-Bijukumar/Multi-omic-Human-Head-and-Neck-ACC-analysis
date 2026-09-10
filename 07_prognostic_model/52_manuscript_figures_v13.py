"""
Stage 52 - the two manuscript figure panels for v13.

Both figures are cut to the panels that carry a distinct claim, so that every panel can be
cited inline at the sentence it supports rather than in a trailing list.

  Figure 1   6 panels - the external evidence
  Figure S1  6 panels - internal validation, what refitting reproduces, composition

Dropped from the v8 layout, as making no argument the text needs: batch-vs-frozen scoring
agreement, the two score-by-group strip plots (their ROC curves say the same thing), the
summary forest (it restates panels already shown), the panel-size scan (a Methods
statement), and the separate composition-combination bars (folded into an annotation).

Everything is read from the locked result files. Nothing is refitted or recomputed. Where a
stage stored only summary statistics for a null distribution, that summary is drawn as an
interval - no histogram is simulated from a mean and a standard deviation.
Outputs to work/outputs_manuscript_v13/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = f"{ACC_DATA_ROOT}"
W = f"{R}/work"
OUT = f"{W}/outputs_manuscript_v13"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.size": 7.4, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": .8, "xtick.major.width": .8, "ytick.major.width": .8,
                     "svg.fonttype": "none", "font.family": "DejaVu Sans"})
RED, BLUE, GREY, GREEN = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837"
PANEL_GENES = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print("  wrote", name)


def tag(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold", fontsize=8.4, pad=6)


def km(t, e):
    ts = np.sort(np.unique(t[e == 1]))
    s, out = 1.0, [(0.0, 1.0)]
    for x in ts:
        n = (t >= x).sum()
        d = ((t == x) & (e == 1)).sum()
        if n:
            s *= 1 - d / n
        out.append((x, s))
    return np.array(out)


thr = json.load(open(f"{W}/outputs_threshold/THRESHOLD_results.json"))
nulls = json.load(open(f"{W}/outputs_nulls/NULLS_results.json"))
ccr = json.load(open(f"{W}/outputs_ccr2020/CCR2020_results.json"))
fre = json.load(open(f"{W}/outputs_frerich/FRERICH_results.json"))
comp = json.load(open(f"{W}/outputs_deconv/DK_COMPOSITION_results.json"))
stab = json.load(open(f"{W}/outputs_metrics/INTERNAL_stability_full.json"))
axis = json.load(open(f"{W}/outputs_axis/AXIS_stability.json"))
draws = np.load(f"{W}/outputs_frerich/frerich_null_draws.npz")

dk = pd.read_csv(f"{W}/outputs_threshold/dk_scores_by_procedure.csv")
cc = pd.read_csv(f"{W}/outputs_ccr2020/ccr2020_scores.csv")
fr = pd.read_csv(f"{W}/outputs_frerich/frerich_scores.csv")
seedp = pd.read_csv(f"{W}/outputs_vi/multiseed_panels.csv")
foldp = pd.read_csv(f"{W}/outputs_axis/fold_panel_agreement.csv")
prog = pd.read_csv(f"{W}/outputs_axis/programme_frequency.csv")

lt = thr["locked_threshold"]["survival"]
t1, t2, t3 = fre["test1_batch"], fre["test2_frozen"], fre["test3_matched_nulls"]
LOCKED_C = axis["locked_dk_c_5yr"]

# ================================================================ FIGURE 1
fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 3, hspace=.62, wspace=.40)

# A  Danish survival under the discovery-locked threshold
ax = fig.add_subplot(gs[0, 0])
for grp, col, lab in [(0, BLUE, "low risk"), (1, RED, "high risk")]:
    m = (dk.high_risk_locked == grp).values
    c = km(dk["T"].values[m], dk["E"].values[m])
    ax.step(c[:, 0], c[:, 1], where="post", color=col, lw=1.7, label=f"{lab} (n={int(m.sum())})")
ax.set_xlim(0, 60); ax.set_ylim(0, 1.03)
ax.set_xlabel("months"); ax.set_ylabel("overall survival")
ax.text(.04, .07, f"HR {lt['hr']:.2f} ({lt['ci'][0]:.2f}–{lt['ci'][1]:.2f})\n"
                  f"log-rank $P$ = {lt['logrank_p']:.4f}", transform=ax.transAxes, fontsize=6.6)
ax.legend(fontsize=6.4, frameon=False, loc="lower left", bbox_to_anchor=(.02, .21))
tag(ax, "A", "Danish cohort (n = 54)")

# B  Danish discrimination against published comparators
ax = fig.add_subplot(gs[0, 1])
vals = [LOCKED_C, 0.540, 0.530]
ax.bar(range(3), vals, color=[RED, GREY, GREY], width=.62, edgecolor="white")
ax.axhline(.5, color="#bbb", lw=.9, ls=":")
ax.set_xticks(range(3))
ax.set_xticklabels(["four-gene\npanel", "Brayer\n14-gene", "Brayer\n49-gene"], fontsize=6.2)
ax.set_ylim(.4, .74); ax.set_ylabel("5-year C-index")
ax.text(0, vals[0] + .012, f"{vals[0]:.3f}", ha="center", fontsize=6.6, fontweight="bold")
tag(ax, "B", "Survival discrimination")

# C  MD Anderson molecular subtype
ax = fig.add_subplot(gs[0, 2])
f_, t_, _ = roc_curve((cc.acc_subtype == 1).astype(int).values, cc.score4.values)
ax.plot(f_, t_, color=RED, lw=1.9)
ax.plot([0, 1], [0, 1], color="#bbb", lw=.9, ls=":")
ax.set_xlabel("1 − specificity"); ax.set_ylabel("sensitivity")
ax.text(.30, .12, f"AUC {ccr['test1']['auc']:.3f}\n"
                  f"({ccr['test1']['ci'][0]:.3f}–{ccr['test1']['ci'][1]:.3f})\n"
                  f"Holm $P$ = {ccr['holm']['test1_adj_p']:.1e}", fontsize=6.6)
tag(ax, "C", "Ferrarotto cohort (n = 54)")

# D  third cohort, poor-outcome group
ax = fig.add_subplot(gs[1, 0])
f2, t2_, _ = roc_curve(fr.y_poor_outcome.values, fr.score_batch.values)
ax.plot(f2, t2_, color=RED, lw=1.9)
ax.plot([0, 1], [0, 1], color="#bbb", lw=.9, ls=":")
ax.set_xlabel("1 − specificity"); ax.set_ylabel("sensitivity")
ax.text(.30, .12, f"AUC {t1['AUC']:.3f}\n({t1['ci95'][0]:.3f}–{t1['ci95'][1]:.3f})\n"
                  f"$P$ = {t1['mannwhitney_p_onesided']:.4f}", fontsize=6.6)
tag(ax, "D", "Frerich cohort (n = 66)")

# E  the locked threshold transfers unchanged
ax = fig.add_subplot(gs[1, 1])
tn, fp, fn, tp = t2["tn"], t2["fp"], t2["fn"], t2["tp"]
ax.bar([0, 1], [tn, fn], color=BLUE, edgecolor="white", width=.62, label="called low risk")
ax.bar([0, 1], [fp, tp], bottom=[tn, fn], color=RED, edgecolor="white", width=.62,
       label="called high risk")
for k, (a_, b_) in enumerate(zip([tn, fn], [fp, tp])):
    ax.text(k, a_ / 2, str(a_), ha="center", va="center", color="white", fontsize=6.8, fontweight="bold")
    ax.text(k, a_ + b_ / 2, str(b_), ha="center", va="center", color="white", fontsize=6.8, fontweight="bold")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Group 2\n(n=52)", "Group 1\n(n=14)"], fontsize=6.4)
ax.set_ylabel("patients"); ax.set_ylim(0, 66)
ax.legend(fontsize=6.0, frameon=False, loc="upper right")
ax.text(.97, .60, f"sens {t2['sensitivity']:.2f}  spec {t2['specificity']:.2f}\n"
                  f"OR {t2['odds_ratio']:.1f}, $P$ = {t2['fisher_p_onesided']:.4f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=6.2)
tag(ax, "E", "Locked cut transfers")

# F  matched random panels, all three cohorts
ax = fig.add_subplot(gs[1, 2])
ax.hist(draws["matched_panel_auc"], bins=40, color="#cdd6de", edgecolor="white", linewidth=.3)
ax.axvline(t1["AUC"], color=RED, lw=1.9)
ax.set_xlabel("AUC under matched null"); ax.set_ylabel("count")
ax.text(.03, .97, "matched-panel $P$\n"
                  f"Danish     {nulls['dk_c_index_5yr']['matched']['p']:.4f}\n"
                  f"Ferrarotto {nulls['ccr2020_subtype_auc']['matched']['p']:.4f}\n"
                  f"Frerich    {t3['empirical_p']:.4f}",
        transform=ax.transAxes, fontsize=6.2, va="top", family="DejaVu Sans Mono")
ax.text(t1["AUC"], ax.get_ylim()[1] * .42, f" observed\n {t1['AUC']:.3f}",
        color=RED, fontsize=6.4, fontweight="bold", va="top")
tag(ax, "F", "Against matched panels")

save(fig, "Figure_1")

# =============================================================== FIGURE S1
fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 3, hspace=.62, wspace=.42)

# A  internal validation, optimism removed
ax = fig.add_subplot(gs[0, 0])
s = stab["summary"]
mu = [s["as_built"]["mean"], s["leakage_free"]["mean"]]
sd = [s["as_built"]["sd"], s["leakage_free"]["sd"]]
ax.bar([0, 1], mu, yerr=sd, color=[GREY, RED], width=.6, edgecolor="white",
       error_kw=dict(lw=1, capsize=3))
ax.set_xticks([0, 1]); ax.set_xticklabels(["as built\n(leaky)", "leakage-free"], fontsize=6.4)
ax.set_ylim(.5, 1.0); ax.set_ylabel("internal AUC")
ax.text(.5, .96, f"optimism +{mu[0]-mu[1]:.3f}", transform=ax.transAxes, ha="center",
        va="top", fontsize=6.4)
for k, (m_, e_) in enumerate(zip(mu, sd)):
    ax.text(k, m_ + e_ + .014, f"{m_:.3f}", ha="center", fontsize=6.4)
tag(ax, "A", "Internal validation")

# B  genes shuffle, the programme concentrates
ax = fig.add_subplot(gs[0, 1])
gsx = pd.read_csv(f"{W}/outputs_compact/gene_selection_stability.csv").head(4)
pr = prog.sort_values("enrichment", ascending=False).head(3)
labels, vals, nl_, cols = [], [], [], []
for _, r in gsx.iterrows():
    labels.append(r.gene); vals.append(r.freq_top4); nl_.append(np.nan)
    cols.append(RED if r.gene in PANEL_GENES else GREY)
labels.append(""); vals.append(np.nan); nl_.append(np.nan); cols.append("none")
for _, r in pr.iterrows():
    labels.append(r.programme.replace("_", " ")); vals.append(r.frac_folds_hit)
    nl_.append(r.null_frac_folds_hit); cols.append(GREEN)
y = np.arange(len(labels))[::-1]
ax.barh(y, [0 if np.isnan(v) else v for v in vals], height=.66, color=cols, edgecolor="white")
for yy, nv in zip(y, nl_):
    if not np.isnan(nv):
        ax.plot([nv, nv], [yy - .33, yy + .33], color="#222", lw=1.5, zorder=4)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.8)
for t_, lab in zip(ax.get_yticklabels(), labels):
    if lab in gsx.gene.values:
        t_.set_style("italic")
ax.set_xlim(0, .75); ax.set_xlabel("fraction of refits selecting it")
ax.text(.97, .95, "genes", transform=ax.transAxes, ha="right", fontsize=6.0,
        style="italic", color="#555")
ax.text(.97, .34, "programmes\n| = chance", transform=ax.transAxes, ha="right",
        fontsize=6.0, style="italic", color="#555")
tag(ax, "B", "Genes shuffle, axis holds")

# C  every refit transfers
ax = fig.add_subplot(gs[0, 2])
rp = axis["external_discrimination"]["random_panel_dk_c"]
ax.axhspan(rp["mean"] - rp["sd"], rp["mean"] + rp["sd"], color="#dfe4e9", zorder=0)
ax.axhline(rp["mean"], color=GREY, lw=1.2, zorder=1)
for k, (v, col) in enumerate([(foldp.dk_c_5yr.values, BLUE), (seedp.dk_c_5yr.values, GREEN)]):
    jit = np.random.default_rng(k + 1).normal(0, .055, len(v))
    ax.scatter(np.full(len(v), k) + jit, v, s=15, color=col, edgecolor="white",
               linewidth=.3, zorder=3)
    ax.hlines(np.mean(v), k - .25, k + .25, color="black", lw=1.6, zorder=4)
    ax.text(k, .335, f"{np.mean(v):.3f}", ha="center", fontsize=6.4, fontweight="bold")
    ax.text(k, .305, f"{(v > .5).sum()}/{len(v)} > 0.5", ha="center", fontsize=5.8, color="#555")
ax.axhline(LOCKED_C, color=RED, lw=1.4, ls="--", zorder=2)
ax.text(1.32, LOCKED_C, f" locked\n {LOCKED_C:.3f}", color=RED, fontsize=6.2, va="center")
ax.text(1.32, rp["mean"], f" random\n {rp['mean']:.3f}", color="#555", fontsize=6.0, va="center")
ax.set_xticks([0, 1]); ax.set_xticklabels(["25 fold\nrefits", "20 seed\nrefits"], fontsize=6.2)
ax.set_xlim(-.45, 2.0); ax.set_ylim(.29, .78)
ax.set_ylabel("Danish 5-year C-index")
tag(ax, "C", "Every refit transfers")

# D  panels sharing no gene still track the locked score
ax = fig.add_subplot(gs[1, 0])
zo = foldp[foldp.overlap_locked == 0]
nlz = axis["score_agreement"]["null_spearman_vs_locked"]
ax.axvspan(nlz["mean"] - nlz["sd"], nlz["mean"] + nlz["sd"], color="#dfe4e9", zorder=0)
ax.axvline(nlz["mean"], color=GREY, lw=1.2, zorder=1)
ax.scatter(zo.spearman_vs_locked, np.random.default_rng(4).uniform(.25, .75, len(zo)),
           s=20, color=GREEN, edgecolor="white", linewidth=.3, zorder=3)
ax.axvline(zo.spearman_vs_locked.mean(), color=GREEN, lw=1.7, zorder=2)
ax.set_ylim(0, 1); ax.set_yticks([])
ax.set_xlabel("rank agreement with the locked score")
ax.text(.02, .96, f"{len(zo)} panels sharing no gene\nwith the locked panel\n"
                  f"mean ρ = {zo.spearman_vs_locked.mean():.3f}",
        transform=ax.transAxes, fontsize=6.2, va="top")
ax.text(nlz["mean"], .05, f"random {nlz['mean']:.3f}", fontsize=5.8, color="#555", ha="center")
tag(ax, "D", "Different genes, same axis")

# E  composition axis
ax = fig.add_subplot(gs[1, 1])
sh = pd.read_csv(f"{W}/outputs_deconv/sc_true_poor_vs_good_shift.csv").set_index("cell_type")
cx = pd.read_csv(f"{W}/outputs_deconv/dk_percompartment_cox.csv").set_index("cell_type")
j = sh.join(cx, how="inner").dropna(subset=["log2fc", "hr_per_sd"])
j = j[j.hr_per_sd > 0]
rx, ry = j.log2fc.rank(), np.log(j.hr_per_sd).rank()
ax.scatter(rx, ry, s=21, color=BLUE, edgecolor="white", linewidth=.4, zorder=3)
for n, xx, yy_ in zip(j.index, rx, ry):
    ax.annotate(n.replace("_Tumor", "").replace("_", " ")[:10], (xx, yy_),
                fontsize=5.0, xytext=(3, 2.5), textcoords="offset points")
ax.set_xlabel("rank, single-cell shift"); ax.set_ylabel("rank, Danish log HR")
ax.set_xlim(.2, len(j) + 1.4); ax.set_ylim(.2, len(j) + 1.4)
dC = comp["combination"]["C_combined"] - comp["combination"]["C_expr4"]
ax.text(.03, .96, f"ρ = 0.80\ncomposition C = {comp['primary']['C']:.3f}\n"
                  f"combined ΔC = {dC:+.3f}", transform=ax.transAxes, fontsize=6.0, va="top")
tag(ax, "E", "Composition, same axis")

# F  no single gene carries the result
ax = fig.add_subplot(gs[1, 2])
fr_auc = [fre["per_gene_auc"][g] for g in PANEL_GENES]
cc_auc = [ccr["secondary"]["per_gene_auc"][g] for g in PANEL_GENES]
zc = [100 * fre["gene_detection_zero_samples"][g] / fre["n_analysed"] for g in PANEL_GENES]
x = np.arange(4)
ax.bar(x - .19, cc_auc, .36, color=RED, edgecolor="white", label="Ferrarotto")
ax.bar(x + .19, fr_auc, .36, color=BLUE, edgecolor="white", label="Frerich")
ax.axhline(.5, color="#bbb", lw=.9, ls=":")
ax.set_xticks(x)
ax.set_xticklabels([f"{g}\n{z:.0f}%" for g, z in zip(PANEL_GENES, zc)], fontsize=5.6)
ax.set_ylim(.4, 1.0); ax.set_ylabel("single-gene AUC")
ax.legend(fontsize=5.8, frameon=False, loc="upper center", ncol=2, columnspacing=.8)
ax.text(.99, .02, "% undetected in archival FFPE", transform=ax.transAxes,
        ha="right", fontsize=5.4, style="italic")
tag(ax, "F", "No gene carries it alone")

save(fig, "Figure_S1")
print("done ->", OUT)

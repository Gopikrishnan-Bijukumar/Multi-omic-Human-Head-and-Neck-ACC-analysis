"""STAGE 7d - deconvolution figures for the paper.

The single most conspicuous gap in the figure set: nothing anywhere plots the scCODA
compositional result against the deconvolution-estimated shift, even though
`work/deconvolv_plan.md` names scCODA as the source of the compositional prior. This script
fills that gap and adds the composition detail figures that were never drawn.

Figures
  M4  deconvolution concordance
        (A) scCODA credible effects - discovery cohort, TRUE counted cells, 3 credible of 26
        (B) true single-cell log2FC vs internal-cohort DWLS-estimated log2FC, 9 compartments
        (C) true single-cell log2FC vs DK per-compartment log2 HR - does the single-cell
            finding transfer to an external bulk cohort?
  S2  composition detail
        (A) per-DK-patient stacked composition, ordered by composition score
        (B) the 3 key compartments by DK 5-year outcome
        (C) internal cohort Poor vs Good, true beside estimated
  S3  gate transparency
        (A) G1/G2/G3 for all three methods against their pre-declared thresholds
        (B) per-patient composition r, making the FAILED G3 (0.685 vs 0.70) visible

FRAMING, carried from section 7 of the consolidated document and not to be softened:
the composition test PASSED but its premise was REFUTED. Composition is prognostic in DK
(HR/SD 1.63, p = 0.039) yet correlates rho = +0.477 with the 4-gene score and adds nothing
(delta-C = -0.006). It is convergent validation, not a second axis. No individual compartment
survives FDR in DK (best q = 0.298). Gate G3 failed, so only rank-based claims are supported.

PROVENANCE NOTE: the canonical scCODA result (annot_5, 24 patients, 26 populations) exists on
disk only as prose, SVG images and executed notebook outputs - the notebook ran with
`csv exports: False`. The only scCODA CSVs on disk are annot_6, a DIFFERENT, non-canonical run
(19 cell types, only `ADC - Tumor` credible) that disagrees with the paper record. This script
parses the stored notebook outputs into sccoda_annot5_results.csv so the figure is driven by a
machine-readable artifact for the first time. The stale annot_6 CSVs are left untouched.

Terminology: the discovery cohort is the "internal cohort" in every label.

Reads only; writes to work/outputs_metrics/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json
import os
import re

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

ROOT = f"{ACC_DATA_ROOT}"
PROJ = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_metrics"
DECONV = f"{ACC_DATA_ROOT}/model_outputs_deconv"
NB = (f"{PROJ}/scripts/named_scripts/sccoda_analysis/"
      "scCODA_annot_5_clinical_outcome_final.ipynb")

# Variant A palette - these figures sit beside GATE_* and DK_composition_*
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "none"})
RED, BLUE, GREY, GREEN = "#b2182b", "#2166ac", "#9a9a9a", "#1b7837"
PURPLE = "#762a83"


def save(fig, name):
    """Every figure in both formats; PNG at 600 dpi, SVG with live text."""
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def panel_label(ax, letter):
    ax.text(-0.16, 1.07, letter, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="left")


# ================================================================ scCODA annot_5 extraction
def parse_sccoda_notebook(path):
    """Recover the annot_5 results table from the notebook's stored outputs.

    The output is a wrapped pandas repr in three column-blocks; every data line begins with the
    row index, so the blocks are re-joined on that.
    """
    nb = json.load(open(path))
    txt = ""
    for o in nb["cells"][10].get("outputs", []):
        if "text" in o:
            txt += "".join(o["text"])
        elif "data" in o and "text/plain" in o["data"]:
            txt += "".join(o["data"]["text/plain"])

    rows = {}
    for ln in txt.splitlines():
        m = re.match(r"^\s*(\d+)\s+(.*\S)\s*$", ln)
        if not m:
            continue
        i, rest = int(m.group(1)), m.group(2)
        rows.setdefault(i, []).append(rest)

    recs = []
    for i in sorted(rows):
        parts = rows[i]
        if len(parts) != 3:
            raise RuntimeError(f"scCODA parse: row {i} has {len(parts)} blocks, expected 3")
        # block 1: "<Cell Type>   <inclusion_probability>"
        b1 = parts[0].rsplit(None, 1)
        ct, incl = b1[0].strip(), float(b1[1])
        # block 2: log2fc, fold_change, hdi_3, hdi_97
        b2 = [float(v) for v in parts[1].split()]
        # block 3: posterior_sd, credible, direction (direction has spaces)
        b3 = parts[2].split(None, 2)
        recs.append({"cell_type": ct, "inclusion_probability": incl,
                     "log2fc": b2[0], "fold_change_vs_Good": b2[1],
                     "hdi_3": b2[2], "hdi_97": b2[3],
                     "posterior_sd": float(b3[0]), "credible": b3[1] == "True",
                     "direction": b3[2].strip()})
    df = pd.DataFrame(recs)
    assert len(df) == 26, f"expected 26 populations, parsed {len(df)}"
    assert df.credible.sum() == 3, f"expected 3 credible, parsed {int(df.credible.sum())}"
    return df


sc5 = parse_sccoda_notebook(NB)
sc5["source"] = "scCODA annot_5, 24 patients, 95441 cells, 100000 draws, seed 20260729"
sc5.to_csv(f"{OUT}/sccoda_annot5_results.csv", index=False)
print("scCODA annot_5 recovered from notebook outputs:")
for _, r in sc5[sc5.credible].iterrows():
    print(f"  {r.cell_type:<32s} log2FC {r.log2fc:+.3f}  incl {r.inclusion_probability:.4f}  "
          f"{r.direction}")
print(f"  -> {OUT}/sccoda_annot5_results.csv")

# ================================================================ inputs
shift = pd.read_csv(f"{DECONV}/sc_true_poor_vs_good_shift.csv")       # true sc, 9 compartments
pcox = pd.read_csv(f"{DECONV}/dk_percompartment_cox.csv")             # DK Cox, 9 compartments
gate = json.load(open(f"{DECONV}/GATE_results.json"))
dk_comp = pd.read_csv(f"{DECONV}/dk_estimated_composition.csv", index_col=0)
jse_est = pd.read_csv(f"{DECONV}/jse_estimated_composition_full.csv", index_col=0)
true_c = pd.read_csv(f"{DECONV}/sc_true_composition_coarse.csv", index_col="patient")
labs = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/sample_labels.csv").set_index("sample")
surv = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs/dk_survival.csv").set_index("UID")

KEY = ["Epithelial_Tumor", "Dividing", "Myoepithelial_Tumor"]
COMPS = list(dk_comp.columns)


def nice(c):
    return c.replace("_", " ")


# ================================================================ M4
fig, ax = plt.subplots(1, 3, figsize=(15.6, 4.8))
fig.subplots_adjust(wspace=.42)

# (A) scCODA: inclusion probability for all 26 populations.
# NOTE: the log2FC column is FDR-thresholded - scCODA sets it to ~0 for every non-credible
# population, so plotting it would show 23 empty bars. The inclusion probability is the actual
# evidence and is never thresholded, so that is what is shown; the credible three are annotated
# with their effect size.
a = ax[0]
d = sc5.sort_values("inclusion_probability")
yy = np.arange(len(d))
cols = [(RED if r.log2fc > 0 else BLUE) if r.credible else GREY for _, r in d.iterrows()]
a.barh(yy, d.inclusion_probability.values, color=cols, height=.74)
a.axvline(0.928, color="k", ls="--", lw=1.1)
a.text(0.915, 17.0, "FDR 5% threshold 0.928", fontsize=6.6, ha="center", va="center",
       rotation=90)
a.axvspan(0.37, 0.54, color=GREY, alpha=.13, zorder=0)
a.text(0.455, 12.0, "no-information range", fontsize=6.4, ha="center", color="#666",
       rotation=90)
for i, (_, r) in enumerate(d.iterrows()):
    if r.credible:
        a.text(r.inclusion_probability - .02, i, f"log2FC {r.log2fc:+.2f}", va="center",
               ha="right", fontsize=6.8, color="white", fontweight="bold")
a.set_yticks(yy)
a.set_yticklabels([c if len(c) < 32 else c[:30] + "..." for c in d.cell_type], fontsize=6.3)
a.set_xlim(0, 1.06)
a.set_ylim(-1.0, len(d) - .3)
a.set_xlabel("posterior inclusion probability")
a.set_title("scCODA, internal cohort\n24 patients, TRUE counted cells - 3 credible of 26",
            fontsize=8.6)
a.legend(handles=[Line2D([], [], color=RED, lw=6, label="credible, enriched in Poor"),
                  Line2D([], [], color=BLUE, lw=6, label="credible, enriched in Good"),
                  Line2D([], [], color=GREY, lw=6, label="not credible")],
         fontsize=6.6, frameon=False, loc="lower right")
panel_label(a, "A")

# (B) true single-cell vs internal-cohort deconvolution-estimated shift
b = ax[1]
jl = labs.loc[[p for p in jse_est.index if p in labs.index], "outcome"]
je = jse_est.loc[jl.index]
est_lfc = {}
for c in COMPS:
    pv = je.loc[jl == "Poor", c].mean()
    gv = je.loc[jl == "Good", c].mean()
    est_lfc[c] = float(np.log2((pv + 1e-4) / (gv + 1e-4)))
sh = shift.set_index("cell_type")
xs = np.array([sh.loc[c, "log2fc"] for c in COMPS])
ys = np.array([est_lfc[c] for c in COMPS])
off = {"Epithelial_Tumor": (8, -12), "Dividing": (8, 4), "Myoepithelial_Tumor": (8, 4)}
for c, x, yv in zip(COMPS, xs, ys):
    k = c in KEY
    b.plot(x, yv, "o", ms=11 if k else 6.5, color=(RED if x > 0 else BLUE) if k else GREY,
           zorder=3 if k else 2)
    if k:
        b.annotate(nice(c), (x, yv), fontsize=7.4, xytext=off[c], textcoords="offset points")
lim = [min(xs.min(), ys.min()) - .5, max(xs.max(), ys.max()) + .5]
b.plot(lim, lim, ls=":", color=GREY, lw=1)
b.axhline(0, color="k", lw=.6)
b.axvline(0, color="k", lw=.6)
b.set_xlim(*lim)
b.set_ylim(*lim)
rb = stats.spearmanr(xs, ys)
b.set_xlabel("true single-cell log2FC (Poor vs Good)")
b.set_ylabel("deconvolution-estimated log2FC")
b.set_title(f"Does bulk deconvolution recover the shift?\n"
            f"internal cohort, 9 compartments\nrho = {rb.correlation:+.3f}", fontsize=8.6)
panel_label(b, "B")

# (C) true single-cell shift vs DK per-compartment hazard.
# Muscle is a 2%-of-cells compartment whose Cox fit is degenerate (HR = 0.006, log2 = -7.4).
# It is kept in the rank correlation - Spearman is rank-based, so it does not distort rho - but
# drawn at the axis edge so the other eight compartments remain readable.
c_ = ax[2]
pc = pcox.set_index("cell_type")
xs2 = np.array([sh.loc[c, "log2fc"] for c in COMPS])
ys2 = np.array([np.log2(pc.loc[c, "hr_per_sd"]) for c in COMPS])
YLO, YHI = -0.95, 0.78
off2 = {"Epithelial_Tumor": (9, 2), "Dividing": (9, -13), "Myoepithelial_Tumor": (9, 3)}
for c, x, yv in zip(COMPS, xs2, ys2):
    k = c in KEY
    clipped = yv < YLO
    yp = YLO + .06 if clipped else yv
    c_.plot(x, yp, "v" if clipped else "o", ms=11 if k else 6.5,
            color=(RED if x > 0 else BLUE) if k else GREY, zorder=3 if k else 2)
    if clipped:
        c_.annotate(f"{nice(c)}  log2HR = {yv:.1f}\n(2% compartment, off scale)", (x, yp),
                    fontsize=6.4, xytext=(9, 2), textcoords="offset points", color="#555")
    elif k:
        c_.annotate(f"{nice(c)}\nq = {pc.loc[c, 'fdr_q']:.2f}", (x, yv), fontsize=7.0,
                    xytext=off2[c], textcoords="offset points")
c_.axhline(0, color="k", lw=.6)
c_.axvline(0, color="k", lw=.6)
c_.set_ylim(YLO, YHI)
c_.set_xlim(xs2.min() - .45, xs2.max() + .75)
rc = stats.spearmanr(xs2, ys2)
keep = np.array([c != "Muscle" for c in COMPS])
rc8 = stats.spearmanr(xs2[keep], ys2[keep])
c_.set_xlabel("true single-cell log2FC (Poor vs Good)")
c_.set_ylabel("DK log2 hazard ratio per SD")
c_.set_title(f"Does it transfer to the external cohort?\n"
             f"DK n = 54, rho = {rc.correlation:+.2f} (all 9) / {rc8.correlation:+.2f} "
             f"(excl. Muscle)", fontsize=8.6)
panel_label(c_, "C")

fig.suptitle("A compositional axis found in single cells, estimated from bulk, and carried to an "
             "external cohort\nall three key compartments move in the predicted direction; "
             "individually none is significant in DK", fontsize=9.6, y=1.07)
save(fig, "M4_deconvolution_concordance")
print("M4 written")

# ================================================================ S2
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.2))

# (A) per-DK-patient stacked composition ordered by composition score
a = ax[0]
zc = (dk_comp - dk_comp.mean(0)) / (dk_comp.std(0) + 1e-12)
TRIAD = {"Epithelial_Tumor": +1.0, "Dividing": +1.0, "Myoepithelial_Tumor": -1.0}
cscore = sum(w * zc[c] for c, w in TRIAD.items())
o = np.argsort(cscore.values)
dd = dk_comp.iloc[o]
order_c = KEY + [c for c in COMPS if c not in KEY]
pal = {"Epithelial_Tumor": RED, "Dividing": "#d6604d", "Myoepithelial_Tumor": BLUE}
base = np.zeros(len(dd))
for c in order_c:
    v = dd[c].values
    a.bar(np.arange(len(dd)), v, bottom=base, width=1.0,
          color=pal.get(c, GREY), alpha=1.0 if c in KEY else .45,
          label=nice(c) if c in KEY else ("other 6 compartments" if c == order_c[3] else None))
    base += v
a.set_xlim(-.5, len(dd) - .5)
a.set_ylim(0, 1)
a.set_xlabel("DK patients, ordered by composition score (low risk -> high risk)")
a.set_ylabel("estimated fraction")
a.set_title("Estimated composition, every DK patient\nn = 54, DWLS", fontsize=8.8)
a.legend(fontsize=6.8, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.17), ncol=4)
panel_label(a, "A")

# (B) key compartments by DK 5-year outcome
b = ax[1]
sv = surv.loc[dk_comp.index]
died5 = ((sv.time.values <= 60) & (sv.event.values == 1))
alive5 = sv.time.values > 60
grp = np.where(died5, "died <5y", np.where(alive5, "alive >5y", "censored"))
pos = 0
for c in KEY:
    for gname, col in [("died <5y", RED), ("alive >5y", BLUE)]:
        v = dk_comp[c].values[grp == gname]
        bp = b.boxplot([v], positions=[pos], widths=.62, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor=col, alpha=.35, edgecolor=col)
        for w in ("whiskers", "caps", "medians"):
            for it in bp[w]:
                it.set(color=col)
        b.plot(np.random.default_rng(0).normal(pos, .07, len(v)), v, "o", ms=3, color=col,
               alpha=.75, zorder=3)
        pos += 1
    pos += .8
b.set_xticks([0.5, 3.3, 6.1])
b.set_xticklabels([nice(c) for c in KEY], fontsize=7.6)
b.set_ylabel("estimated fraction")
b.set_title("DK compartments by 5-year outcome\ndirection as predicted; none survives FDR",
            fontsize=8.8)
b.legend(handles=[Line2D([], [], color=RED, lw=6, alpha=.5, label="died <5y"),
                  Line2D([], [], color=BLUE, lw=6, alpha=.5, label="alive >5y")],
         fontsize=7.2, frameon=False, loc="upper right")
panel_label(b, "B")

# (C) internal cohort, true beside estimated
c_ = ax[2]
tc = true_c.loc[[p for p in true_c.index if p in labs.index]]
tl = labs.loc[tc.index, "outcome"]
w = .34
xx = np.arange(len(KEY))
for i, (src, frame, lab_, alpha) in enumerate(
        [("true", tc, "true (counted cells)", 1.0), ("est", je, "estimated (bulk)", .55)]):
    ll = tl if src == "true" else jl
    poor = [frame.loc[ll == "Poor", c].mean() for c in KEY]
    good = [frame.loc[ll == "Good", c].mean() for c in KEY]
    c_.bar(xx + (i - .5) * w - .09, poor, w * .46, color=RED, alpha=alpha,
           label=f"Poor, {lab_}")
    c_.bar(xx + (i - .5) * w + .09, good, w * .46, color=BLUE, alpha=alpha,
           label=f"Good, {lab_}")
c_.set_xticks(xx)
c_.set_xticklabels([nice(c) for c in KEY], fontsize=7.6)
c_.set_ylabel("mean fraction")
c_.set_title("Internal cohort: true vs estimated\nsame direction, compressed magnitude "
             "(G3 failed)", fontsize=8.8)
c_.legend(fontsize=6.6, frameon=False, ncol=1)
panel_label(c_, "C")

save(fig, "S2_composition_detail")
print("S2 written")

# ================================================================ S3
fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0))

a = ax[0]
methods = ["nusvr", "dwls", "nnls"]
g = gate["gate"]
xx = np.arange(len(methods))
w = .26
g1 = [g[m]["median_spearman"] for m in methods]
g3 = [g[m]["mean_patient_r"] for m in methods]
g2 = [min(g[m]["dominant"][c] for c in KEY) for m in methods]
a.bar(xx - w, g1, w, color=BLUE, label="G1  median per-cell-type rho")
a.bar(xx, g2, w, color=GREEN, label="G2  worst of the 3 key compartments")
a.bar(xx + w, g3, w, color=PURPLE, label="G3  mean per-patient r")
a.axhline(.40, color=GREY, ls="--", lw=1)
a.text(-0.42, .415, "G1/G2 threshold 0.40", fontsize=6.8, color=GREY, va="bottom")
a.axhline(.70, color=RED, ls="--", lw=1)
a.text(-0.42, .715, "G3 threshold 0.70", fontsize=6.8, color=RED, va="bottom")
for i, m in enumerate(methods):
    if np.isfinite(g3[i]):
        a.text(i + w, g3[i] + .012, f"{g3[i]:.3f}", ha="center", fontsize=7,
               color=RED if g3[i] < .70 else "k")
a.set_xlim(-.55, len(methods) - .45)
a.set_xticks(xx)
a.set_xticklabels([m.upper() for m in methods])
a.set_ylim(0, .95)
a.set_ylabel("gate statistic")
a.set_title("Pre-declared gates, all three methods\nDWLS wins; G3 fails for every method",
            fontsize=8.8)
a.legend(fontsize=7, frameon=False, loc="upper left")
panel_label(a, "A")

b = ax[1]
best = gate["best_method"]
est = pd.read_csv(f"{DECONV}/jse_estimated_composition_{best}.csv", index_col=0)
common = [p for p in est.index if p in true_c.index]
rs = [float(np.corrcoef(est.loc[p, COMPS].values.astype(float),
                        true_c.loc[p, COMPS].values.astype(float))[0, 1]) for p in common]
b.hist(rs, bins=12, color=PURPLE, alpha=.55, edgecolor="white")
b.axvline(np.mean(rs), color=PURPLE, lw=2, label=f"mean = {np.mean(rs):.3f}")
b.axvline(.70, color=RED, ls="--", lw=1.5, label="G3 threshold = 0.70")
b.set_xlabel("per-patient composition correlation (estimated vs true)")
b.set_ylabel("patients")
b.set_title(f"The gate that failed\n{best.upper()}, {len(common)} internal-cohort patients, "
            f"leave-one-out", fontsize=8.8)
b.legend(fontsize=7.4, frameon=False)
panel_label(b, "B")

fig.suptitle("Deconvolution was accepted with a declared deviation, not a passed gate",
             fontsize=9.6, y=1.03)
save(fig, "S3_gate_transparency")
print("S3 written")

print(f"\nfigures written -> {OUT}")

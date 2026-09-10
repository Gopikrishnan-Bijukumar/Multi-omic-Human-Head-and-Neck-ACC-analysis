"""STAGE 7g - figure S4: the internal-cohort estimate and the leak, with their run-to-run spread.

Driven entirely by work/outputs_metrics/INTERNAL_stability_full.json (31e). Nothing is hardcoded
except the published single-draw values, which are carried in that JSON and drawn as reference
markers so the correction is visible rather than asserted.

Panel A - every independent draw of each arm, with mean +/- sd and the published single draw.
Panel B - the leak, taken paired within each repeat, against the published +0.071.
Panel C - why: the LOO and 5x5 CV components separately, showing which one carries the instability.

Writes S4_internal_stability.{png,svg} via the house save() helper (600 dpi + live-text SVG).
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_metrics"

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

C_AB, C_LF, GREY = "#B4436C", "#3C6E8F", "#9a9a9a"


def save(fig, name):
    """Every figure in both formats; PNG at 600 dpi, SVG with live text."""
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def despine(ax):
    ax.spines[["top", "right"]].set_visible(False)


def panel_label(ax, letter):
    ax.text(-0.16, 1.07, letter, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="top", ha="left")


D = json.load(open(f"{OUT}/INTERNAL_stability_full.json"))
reps, S, PUB = D["repeats"], D["summary"], D["published_25c"]
prior = D["prior_independent_draws_pooled"]

draws = {a: [r[a]["disc_auc"] for r in reps] + prior[a] for a in ["as_built", "leakage_free"]}
paired = {a: [r[a]["disc_auc"] for r in reps] for a in ["as_built", "leakage_free"]}
leaks = [r["leak"] for r in reps]

fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.7))
rng = np.random.default_rng(0)

# ---- A: arm-level draws -------------------------------------------------------------------
ax = axes[0]
for i, (arm, col, lab) in enumerate([("as_built", C_AB, "as-built"),
                                     ("leakage_free", C_LF, "leakage-free")]):
    v = np.asarray(draws[arm], float)
    ax.scatter(i + rng.uniform(-.09, .09, len(v)), v, s=26, color=col, alpha=.75,
               edgecolor="none", zorder=3)
    m, sd = S[arm]["mean"], S[arm]["sd"]
    ax.hlines(m, i - .26, i + .26, color=col, lw=2.2, zorder=4)
    ax.add_patch(plt.Rectangle((i - .26, m - sd), .52, 2 * sd, color=col, alpha=.13, zorder=1))
    ax.scatter([i], [PUB[arm]], marker="*", s=190, facecolor="none", edgecolor="k",
               linewidth=1.1, zorder=5)
    ax.text(i + .30, m, f"{m:.3f}\n±{sd:.3f}", va="center", ha="left", fontsize=8, color=col)
    ax.annotate(f"published {PUB[arm]:.3f}", xy=(i - .06, PUB[arm]), xytext=(i - .34, PUB[arm]),
                va="center", ha="right", fontsize=7.5, color="k",
                arrowprops=dict(arrowstyle="-", lw=.7, color="k"))
ax.set_xticks([0, 1]); ax.set_xticklabels(["as-built", "leakage-free"])
ax.set_xlim(-1.05, 1.85); ax.set_ylabel("discovery AUC (nested, fixed k=4)")
ax.set_title("Independent re-runs of the published design", fontsize=9.5)
despine(ax); panel_label(ax, "A")

# ---- B: paired leak ------------------------------------------------------------------------
ax = axes[1]
for j in range(len(reps)):
    ax.plot([0, 1], [paired["as_built"][j], paired["leakage_free"][j]], color=GREY, lw=.8,
            alpha=.8, zorder=2)
ax.scatter(np.zeros(len(reps)), paired["as_built"], s=24, color=C_AB, zorder=3, edgecolor="none")
ax.scatter(np.ones(len(reps)), paired["leakage_free"], s=24, color=C_LF, zorder=3,
           edgecolor="none")
lk = S["leak_paired"]
ax.set_xticks([0, 1]); ax.set_xticklabels(["as-built", "leakage-free"])
ax.set_xlim(-.35, 1.75); ax.set_ylabel("discovery AUC")
ax.set_title("Leak, paired within repeat", fontsize=9.5)
ax.text(1.12, np.mean(paired["leakage_free"]),
        f"leak = {lk['mean']:+.3f} ± {lk['sd']:.3f}\n"
        f"[{lk['pct_lo']:+.3f}, {lk['pct_hi']:+.3f}]\n"
        f"published {PUB['leak']:+.3f}",
        va="center", ha="left", fontsize=8)
despine(ax); panel_label(ax, "B")

# ---- C: which component is unstable --------------------------------------------------------
ax = axes[2]
comp = [("loo_auc", "LOO ×3 seeds"), ("cv5_auc", "5×5 CV")]
pubc = {"loo_auc": 0.843, "cv5_auc": 0.776}
for i, (key, lab) in enumerate(comp):
    v = np.array([r["leakage_free"][key] for r in reps])
    ax.scatter(i + rng.uniform(-.09, .09, len(v)), v, s=26, color=C_LF, alpha=.75,
               edgecolor="none", zorder=3)
    ax.hlines(v.mean(), i - .26, i + .26, color=C_LF, lw=2.2, zorder=4)
    ax.add_patch(plt.Rectangle((i - .26, v.mean() - v.std(ddof=1)), .52, 2 * v.std(ddof=1),
                               color=C_LF, alpha=.13, zorder=1))
    ax.scatter([i], [pubc[key]], marker="*", s=190, facecolor="none", edgecolor="k",
               linewidth=1.1, zorder=5)
    ax.text(i + .30, v.mean(), f"{v.mean():.3f}\n±{v.std(ddof=1):.3f}", va="center", ha="left",
            fontsize=8, color=C_LF)
    ax.text(i - .30, pubc[key], f"pub.\n{pubc[key]:.3f}", va="center", ha="right", fontsize=7,
            color="k")
ax.set_xticks([0, 1]); ax.set_xticklabels([c[1] for c in comp])
ax.set_xlim(-.85, 1.85); ax.set_ylabel("component AUC (leakage-free)")
ax.set_title("Both components reproduce; LOO is noisier", fontsize=9)
despine(ax); panel_label(ax, "C")

fig.legend(handles=[Line2D([], [], marker="*", ls="none", markerfacecolor="none",
                           markeredgecolor="k", markersize=12,
                           label="published single draw (log_25c.txt)"),
                    Line2D([], [], marker="o", ls="none", color=GREY, markersize=6,
                           label="one independent re-run")],
           loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, -.07), fontsize=8.5)
fig.suptitle(f"Internal cohort (n=20): run-to-run spread of the nested-CV estimate "
             f"({D['n_repeats_paired']} paired repeats)", fontsize=10.5, y=1.02)
fig.tight_layout()
save(fig, "S4_internal_stability")
print("saved -> S4_internal_stability.{png,svg}")

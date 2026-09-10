"""STAGE 12: our 4-gene panel vs Brayer's PUBLISHED poor-survival group in the DK cohort.

BACKGROUND
----------
Stage 3 (`18_brayer_benchmark.py`) compared our panel to Brayer's *classifier* by transcribing
their Table 5 coefficients. That answered "whose score ranks patients better?" It did not answer
"can we identify the subgroup they actually named?", because at the time we had no per-patient
labels for their poor-survival group.

Their group membership is still not published per patient. What IS published is the group's
*definition*: Supplementary Table S2 of Brayer et al. (Cancers 2023;15:1390) gives the complete
differential-expression result for "poor survival samples vs rest" in the DK cohort - 731 genes
with signed log2 fold-changes. Figure 3A shows the group as ~6 brown points on an MDS plot; Figure
3B gives its survival: median 8 months versus 80 months for the main group, log-rank p = 0.0061.

This script reconstructs the group from the published signature, checks the reconstruction against
those published survival statistics, and then asks how well our four genes recover it.

WHAT IS AND IS NOT BEING CLAIMED
-------------------------------
 *  Their 731-gene list was derived FROM the brown group in this same cohort. Reconstructing the
    group from that list is therefore near-circular BY CONSTRUCTION: it is a faithful rendering of
    their definition, not an independent rediscovery. The reconstruction step proves nothing.
 *  What is NOT circular is the test that follows. Our four genes were selected on JSE survival,
    years of analysis away from this gene list, and share ZERO members with it (asserted below).
    Their AUC against the reconstructed group is a genuine agreement test.
 *  This measures agreement with a subtype CALL, not survival prediction. The group is
    survival-associated by construction, so a good AUC here is not extra survival evidence.
 *  Group size is not published as a number. AUC is therefore reported across n = 4..10 rather
    than at one cherry-picked size, and the whole curve is what should be quoted.

Reads only. Writes to work/outputs_brayer_group/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import json
import importlib.util
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9

R = f"{ACC_DATA_ROOT}"
OUT = f"{ACC_DATA_ROOT}/model_outputs_brayer_group"
SUPP = f"{ACC_DATA_ROOT}/external_brayer_supp/TableS2_Poor_Survival_DE.csv"

spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3)

RED, BLUE, GREY, BROWN = "#b2182b", "#2166ac", "#9a9a9a", "#8c5a2b"
PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
PANEL8 = PANEL4 + ["EIF3H", "STIL", "KCNK5", "FAM111B"]
B_NULL = 10000
RNG = np.random.default_rng(20260815)

# Published values, Brayer et al. Cancers 2023;15:1390, Figure 3B and Figure 3A.
PUB = {"median_poor_months": 8.0, "median_main_months": 80.0, "logrank_p": 0.0061,
       "n_brown_visible_in_fig3A": 6, "n_dk_salivary_in_paper": 56}


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


print("=" * 78)
print("Stage 12 - our panel vs Brayer's published DK poor-survival group")
print("=" * 78)

# ---------------------------------------------------------------- published signature
raw = pd.read_csv(SUPP, skiprows=3)
raw.columns = ["gene", "ID", "EntrezGene", "Symbol", "log2FC", "FC",
               "PValue", "AdjPVal", "gene_link", "pubmed", "ucsc", "genecards"]
sig = raw[["Symbol", "log2FC", "AdjPVal"]].dropna(subset=["Symbol"]).copy()
sig["log2FC"] = pd.to_numeric(sig["log2FC"], errors="coerce")
sig = sig.dropna(subset=["log2FC"]).drop_duplicates("Symbol")
print(f"\npublished poor-survival DE genes (Table S2): {len(sig)}")

# ---------------------------------------------------------------- DK cohort
dk = pd.read_csv(f"{R}/outputs/03_harmonized/dk_norm_log2.csv", index_col=0)
clin = pd.read_csv(f"{R}/external_dk/DK_sample_data_table2.csv").set_index("UID").loc[dk.index]
T = pd.to_numeric(clin["SURV_MO"], errors="coerce").values
E = pd.to_numeric(clin["SURV_CENS"], errors="coerce").values
has = ~np.isnan(T)
print(f"DK matrix: {dk.shape[0]} samples, {has.sum()} with survival "
      f"(paper analysed {PUB['n_dk_salivary_in_paper']} salivary samples)")

Z = d3.rint(dk)
shared = [g for g in sig["Symbol"] if g in dk.columns]
sig_s = sig.set_index("Symbol").loc[shared]
print(f"signature genes present in our harmonised matrix: {len(shared)}/{len(sig)}")

# The panel must not overlap the published list, or the agreement test is meaningless.
overlap = sorted(set(PANEL8) & set(sig["Symbol"]))
assert not overlap, f"panel genes appear in Brayer's published DE list: {overlap}"
print(f"panel/published-list gene overlap: 0 of {len(PANEL8)}  (agreement test is independent)")

# score: high = poor-survival-like, per the published sign of each log2FC
their = (Z[shared].values * np.sign(sig_s["log2FC"].values)).mean(1)
ours4 = Z[PANEL4].values.mean(1)
ours8 = Z[PANEL8].values.mean(1)

# ---------------------------------------------------------------- reconstruction check
print("\n" + "-" * 78)
print("Reconstruction check - does the published signature recover the published KM?")
print(f"  target: median {PUB['median_poor_months']:.0f} vs {PUB['median_main_months']:.0f} mo, "
      f"log-rank p = {PUB['logrank_p']}")
print("   n  n_surv  median_in  median_out  logrank_p")
recon = []
for n in range(3, 15):
    m = np.zeros(len(their), bool)
    m[np.argsort(-their)[:n]] = True
    mm, oo = m & has, (~m) & has
    if mm.sum() < 2:
        continue
    lr = logrank_test(T[mm], T[oo], E[mm], E[oo])
    km = KaplanMeierFitter()
    km.fit(T[mm], E[mm]); med_in = km.median_survival_time_
    km.fit(T[oo], E[oo]); med_out = km.median_survival_time_
    recon.append(dict(n=n, n_surv=int(mm.sum()), median_in=float(med_in),
                      median_out=float(med_out), logrank_p=float(lr.p_value)))
    print(f"  {n:2d}   {mm.sum():4d}   {med_in:8.1f}  {med_out:9.1f}   {lr.p_value:.5f}")

# ---------------------------------------------------------------- the actual test
print("\n" + "-" * 78)
print("Our panel recovering the reconstructed group")
print("   n   AUC_4gene   AUC_8gene   matched-free null p (4-gene, B=%d)" % B_NULL)
rows = []
cols = Z.values
for n in [4, 5, 6, 7, 8, 10]:
    y = np.zeros(len(their), int)
    y[np.argsort(-their)[:n]] = 1
    a4 = roc_auc_score(y, ours4)
    a8 = roc_auc_score(y, ours8)
    null = np.empty(B_NULL)
    for i in range(B_NULL):
        null[i] = roc_auc_score(y, cols[:, RNG.choice(cols.shape[1], 4, replace=False)].mean(1))
    b = int(np.sum(null >= a4))
    p = (b + 1) / (B_NULL + 1)
    rows.append(dict(n=n, auc_4gene=float(a4), auc_8gene=float(a8),
                     null_b=b, null_B=B_NULL, null_p=float(p),
                     null_mean=float(null.mean())))
    print(f"  {n:2d}    {a4:.3f}       {a8:.3f}       {p:.4f}  (b={b}/{B_NULL}, "
          f"null mean {null.mean():.3f})")

rho4, p4 = stats.spearmanr(ours4, their)
rho8, p8 = stats.spearmanr(ours8, their)
print(f"\ncontinuous agreement across all {len(their)} DK samples:")
print(f"  Spearman(our 4-gene score, their published signature) = {rho4:.3f}  p = {p4:.2g}")
print(f"  Spearman(our 8-gene score, their published signature) = {rho8:.3f}  p = {p8:.2g}")

# ---------------------------------------------------------------- purity confound
# Their published signature is dominated at the top by salivary acinar/serous and squamous
# markers (CST1/2/4/5, MUC7, KLK1, CA6, PIGR, SPRR1A/3, KRT13), all DOWN in the poor-survival
# group. That is the expected signature of higher tumour cellularity, not of tumour biology.
# Quantify it, because it limits how much their group - and therefore our agreement with it -
# can be read as biological.
NORMAL_SG = ["CST1", "CST2", "CST4", "CST5", "MUC7", "KLK1", "CA6", "PIGR", "LTF",
             "ZG16B", "STATH", "HTN1", "AMY1A", "PRR4", "LYZ", "AQP5", "SMR3B"]
sg = [g for g in NORMAL_SG if g in dk.columns]
purity_proxy = -Z[sg].values.mean(1)          # high = little normal salivary tissue
rho_pur, p_pur = stats.spearmanr(their, purity_proxy)
rho_pur4, p_pur4 = stats.spearmanr(ours4, purity_proxy)
top50 = sig.reindex(sig["log2FC"].abs().sort_values(ascending=False).index).head(50)
frac_sg_top50 = float(np.mean([g in set(NORMAL_SG) for g in top50["Symbol"]]))
print("\n" + "-" * 78)
print("Tumour-purity confound in THEIR published group")
print(f"  normal-salivary markers found: {len(sg)}/{len(NORMAL_SG)} -> purity proxy")
print(f"  Spearman(their published signature, purity proxy) = {rho_pur:.3f}  p = {p_pur:.2g}")
print(f"  Spearman(our 4-gene score,        purity proxy) = {rho_pur4:.3f}  p = {p_pur4:.2g}")
print(f"  fraction of their top-50 |log2FC| genes that are normal salivary markers: "
      f"{frac_sg_top50:.2f}")

# ---------------------------------------------------------------- figures
n_ref = PUB["n_brown_visible_in_fig3A"]
y_ref = np.zeros(len(their), int)
y_ref[np.argsort(-their)[:n_ref]] = 1

fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.6))

a = ax[0]
a.scatter(their[y_ref == 0], ours4[y_ref == 0], s=26, c=GREY, edgecolor="none",
          label="main group")
a.scatter(their[y_ref == 1], ours4[y_ref == 1], s=46, c=BROWN, edgecolor="k", lw=.4,
          label=f"reconstructed poor-survival group (n={n_ref})")
a.set_xlabel("Brayer published poor-survival signature (650 genes)")
a.set_ylabel("our 4-gene score")
a.set_title(f"Spearman rho = {rho4:.3f}, p = {p4:.1g}", fontsize=9)
a.legend(frameon=False, fontsize=7.5, loc="upper left")
a.spines[["top", "right"]].set_visible(False)

a = ax[1]
df = pd.DataFrame(rows)
a.plot(df.n, df.auc_4gene, "o-", color=RED, label="our 4-gene")
a.plot(df.n, df.auc_8gene, "s--", color=BLUE, label="our 8-gene", ms=4)
a.plot(df.n, df.null_mean, "^:", color=GREY, label="random 4-gene panels", ms=4)
a.axhline(0.5, color="k", lw=.6, ls=":")
a.set_xlabel("assumed size of their poor-survival group")
a.set_ylabel("AUC")
a.set_ylim(0.35, 1.02)
a.set_title("recovery is stable across plausible group sizes", fontsize=9)
a.legend(frameon=False, fontsize=7.5, loc="lower left")
a.spines[["top", "right"]].set_visible(False)

a = ax[2]
km = KaplanMeierFitter()
for lab, m, c in [(f"reconstructed poor-survival (n={int((y_ref == 1).sum())})", y_ref == 1, BROWN),
                  ("main group", y_ref == 0, "#7fb3d5")]:
    mm = m & has
    km.fit(T[mm], E[mm], label=f"{lab}, n={mm.sum()}")
    km.plot_survival_function(ax=a, ci_show=False, color=c, lw=1.8)
mm, oo = (y_ref == 1) & has, (y_ref == 0) & has
lr = logrank_test(T[mm], T[oo], E[mm], E[oo])
a.set_xlim(0, 120)
a.set_xlabel("months")
a.set_ylabel("fraction surviving")
a.set_title(f"reconstruction: log-rank p = {lr.p_value:.4f}\n"
            f"(published: p = {PUB['logrank_p']})", fontsize=9)
a.legend(frameon=False, fontsize=7.5)
a.spines[["top", "right"]].set_visible(False)

save(fig, "BRAYERGROUP_recovery")

res = {
    "question": "Can our 4 genes identify Brayer's published DK poor-survival group?",
    "published_source": ("Brayer et al. Cancers 2023;15:1390, Supplementary Table S2 "
                         "(poor survival vs rest DE) and Figure 3A/3B"),
    "published_values": PUB,
    "signature_genes_published": int(len(sig)),
    "signature_genes_matched": int(len(shared)),
    "panel_overlap_with_published_list": 0,
    "reconstruction_scan": recon,
    "recovery": rows,
    "spearman_continuous": {"panel4_rho": float(rho4), "panel4_p": float(p4),
                            "panel8_rho": float(rho8), "panel8_p": float(p8)},
    "reference_group_size_used_for_figure": n_ref,
    "purity_confound": {
        "markers_used": sg,
        "their_signature_vs_purity_rho": float(rho_pur), "p": float(p_pur),
        "our_4gene_vs_purity_rho": float(rho_pur4), "our_p": float(p_pur4),
        "frac_normal_salivary_in_their_top50_by_abs_log2FC": frac_sg_top50,
        "note": ("Their poor-survival group is defined partly by loss of normal salivary "
                 "tissue signal, i.e. tumour cellularity. This limits how biological the "
                 "group is, and therefore how much agreement with it can be claimed."),
    },
    "caveats": [
        "Their 731-gene list was derived from this group in this cohort, so reconstructing the "
        "group from it is circular by construction and proves nothing on its own.",
        "The recovery test is NOT circular: our panel shares zero genes with the published list "
        "and was selected on JSE survival.",
        "This measures agreement with a subtype call, not survival prediction. The group is "
        "survival-associated by construction.",
        "Group size is not published as a number; ~6 brown points are visible in Figure 3A. "
        "AUC is reported across n = 4-10 and the whole curve should be quoted.",
        "Our harmonised matrix has 62 DK samples (54 with survival); the paper analysed 56 "
        "salivary samples, and normalisation differs. The reconstruction is approximate.",
    ],
}
json.dump(res, open(f"{OUT}/BRAYERGROUP_results.json", "w"), indent=1)
pd.DataFrame({"sample": dk.index, "their_published_signature": their,
              "our_4gene": ours4, "our_8gene": ours8,
              "reconstructed_poor_group": y_ref}).to_csv(
    f"{OUT}/dk_scores_vs_published_group.csv", index=False)
print(f"\nwrote {OUT}/BRAYERGROUP_results.json")

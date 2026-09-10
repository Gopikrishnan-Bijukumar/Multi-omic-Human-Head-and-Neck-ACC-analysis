"""
Stage 5 — external test of the locked 4-gene panel in the Ferrarotto/Mitani CCR-2020 cohort.

Procedure is fixed by work/PRESPEC_stage5.md, written before this script was run.

TEST 1 (primary)  AUC of the 4-gene score for ACC-I vs ACC-II (20 vs 34).
                  criterion: AUC > 0.5 and Holm-corrected Mann-Whitney p < 0.05.
TEST 2 (primary)  Spearman(4-gene score, TTDeath) among the 34 decedents.
                  criterion: negative rho and Holm-corrected p < 0.05.

Holm is applied across exactly these two primaries. Everything else is secondary and uncorrected.

Nothing is refitted. The panel, its weights and its direction come from outputs/FINAL_model.json.
Reads only; all outputs go to work/outputs_ccr2020/.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"

ROOT = f"{ACC_DATA_ROOT}"
EXT = f"{ACC_DATA_ROOT}/external_ccr2020"
OUT = f"{ACC_DATA_ROOT}/model_outputs_ccr2020"
RNG = np.random.default_rng(20260815)

PANEL4 = ["DSCAM", "ODC1", "NCAPG", "CCNB2"]
PANEL8 = PANEL4 + ["EIF3H", "STIL", "KCNK5", "FAM111B"]


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)


def rint(df):
    """Rank-inverse-normal within cohort, gene-wise across samples (03_data.py:9)."""
    r = df.rank(axis=0)
    n = df.shape[0]
    return pd.DataFrame(stats.norm.ppf((r.values - 0.375) / (n + 0.25)),
                        index=df.index, columns=df.columns)


def boot_auc(y, s, n=2000):
    out = []
    idx = np.arange(len(y))
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        out.append(roc_auc_score(y[b], s[b]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def holm(pvals):
    """Holm-Bonferroni adjusted p-values, order preserved."""
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * pvals[i])
        adj[i] = min(run, 1.0)
    return adj


# ---------------- load ----------------
print("=" * 72)
print("Stage 5 - Ferrarotto/Mitani CCR-2020 external cohort")
print("=" * 72)

expr = pd.read_excel(f"{EXT}/ACC_RNAseq.xlsx", sheet_name="ACC_Mitani_data_RPKM")

# pre-spec step 1: drop Excel-corrupted symbols (dates); step 2: sum RPKM within symbol
keep = expr["gene_name"].map(lambda g: isinstance(g, str))
print(f"dropped Excel date-corrupted gene rows : {int((~keep).sum())}")
expr = expr[keep].copy()
expr["gene_name"] = expr["gene_name"].astype(str)
n_dup = int(expr["gene_name"].duplicated().sum())
expr = expr.groupby("gene_name").sum()
print(f"duplicate symbols collapsed by sum     : {n_dup}")
print(f"expression matrix                      : {expr.shape[0]} genes x {expr.shape[1]} samples")

clin = pd.read_excel(f"{EXT}/CCR2020_Clinical.xlsx", sheet_name="Annotations")
clin = clin.set_index("TID")
samples = [s for s in expr.columns if s in clin.index]
assert len(samples) == 54, len(samples)
expr = expr[samples]
clin = clin.loc[samples]
print(f"matched samples                        : {len(samples)}")

# pre-spec step 3: RINT gene-wise across the 54 samples  (samples on rows)
Z = rint(expr.T)

subtype = clin["ACC Subtype"].astype(int)
y_acc1 = (subtype == 1).astype(int).values          # ACC-I is the positive class
status = clin["Status"].astype(str).values
ttd = pd.to_numeric(clin["TTDeath"], errors="coerce")
dead = status == "D"
print(f"ACC-I / ACC-II                         : {int(y_acc1.sum())} / {int((1 - y_acc1).sum())}")
print(f"decedents with exact death time        : {int(dead.sum())}  "
      f"({np.nanmin(ttd[dead]):.0f}-{np.nanmax(ttd[dead]):.0f} days)")

# ---------------- locked score ----------------
_fm = json.load(open(f"{ACC_DATA_ROOT}/model_outputs/FINAL_model.json"))
beta = dict(zip(_fm["genes"], _fm["beta"]))


def panel_score(genes):
    """mean(sign(beta) * RINT(expr)) -- scheme='sign', 13_compact_panel.py:88."""
    return (Z[genes] * np.sign([beta[g] for g in genes])).mean(axis=1).values


s4 = panel_score(PANEL4)
s8 = panel_score(PANEL8)
print(f"\nlocked panel betas (all positive)      : "
      + ", ".join(f"{g} {beta[g]:+.4f}" for g in PANEL4))

res = {"cohort": "Ferrarotto/Mitani CCR 2020 (MD Anderson, n=54)",
       "source": "Mendeley doi:10.17632/6sbv7bpj5n.1 (CC BY 4.0)",
       "prespec": "work/PRESPEC_stage5.md",
       "panel4": PANEL4, "n_acc1": int(y_acc1.sum()), "n_acc2": int((1 - y_acc1).sum()),
       "n_decedents": int(dead.sum())}

# ================= TEST 1 =================
print("\n" + "=" * 72)
print("TEST 1 (primary)  ACC-I vs ACC-II, locked 4-gene score")
print("=" * 72)

auc1 = roc_auc_score(y_acc1, s4)
mw1 = stats.mannwhitneyu(s4[y_acc1 == 1], s4[y_acc1 == 0], alternative="two-sided")
lo1, hi1 = boot_auc(y_acc1, s4)
print(f"AUC = {auc1:.3f}  [95% CI {lo1:.3f}-{hi1:.3f}]   Mann-Whitney p = {mw1.pvalue:.4g}")

# nulls: random 4-gene panels from genes expressed in this cohort
expressed = expr.index[(expr > 1).mean(axis=1) >= 0.25]
Ze = Z[expressed]
null_auc = []
for _ in range(2000):
    g = RNG.choice(len(expressed), 4, replace=False)
    null_auc.append(roc_auc_score(y_acc1, Ze.values[:, g].mean(1)))
null_auc = np.array(null_auc)
p_rand = float((np.abs(null_auc - 0.5) >= abs(auc1 - 0.5)).mean())

perm = np.array([roc_auc_score(RNG.permutation(y_acc1), s4) for _ in range(2000)])
p_perm = float((np.abs(perm - 0.5) >= abs(auc1 - 0.5)).mean())
print(f"random 4-gene panels (n=2000)  p = {p_rand:.4f}   (null mean AUC {null_auc.mean():.3f})")
print(f"label permutation    (n=2000)  p = {p_perm:.4f}")

res["test1"] = {"auc": float(auc1), "ci": [lo1, hi1], "mw_p": float(mw1.pvalue),
                "random_panel_p": p_rand, "permutation_p": p_perm,
                "n_expressed_universe": int(len(expressed))}

# ================= TEST 2 =================
print("\n" + "=" * 72)
print("TEST 2 (primary)  time to death among decedents, locked 4-gene score")
print("=" * 72)

sd = s4[dead]
td = ttd[dead].values.astype(float)
rho, p_rho = stats.spearmanr(sd, td)
print(f"Spearman rho = {rho:+.3f}   p = {p_rho:.4g}   (n = {len(td)})")
print(f"pre-declared direction is negative -> {'as predicted' if rho < 0 else 'WRONG DIRECTION'}")

early = (td < 5 * 365.25).astype(int)
auc_early = roc_auc_score(early, sd) if len(np.unique(early)) > 1 else np.nan
print(f"secondary: AUC death <5y vs >=5y among decedents = {auc_early:.3f} "
      f"({int(early.sum())} early / {int((1 - early).sum())} late)")

res["test2"] = {"spearman_rho": float(rho), "p": float(p_rho), "n": int(len(td)),
                "direction_as_predicted": bool(rho < 0),
                "secondary_auc_early_death": float(auc_early),
                "n_early": int(early.sum()), "n_late": int((1 - early).sum())}

# ================= Holm across the two primaries =================
adj = holm([mw1.pvalue, p_rho])
res["holm"] = {"test1_adj_p": float(adj[0]), "test2_adj_p": float(adj[1])}
t1_pass = bool(auc1 > 0.5 and adj[0] < 0.05)
t2_pass = bool(rho < 0 and adj[1] < 0.05)
res["test1"]["pass"] = t1_pass
res["test2"]["pass"] = t2_pass

print("\n" + "-" * 72)
print(f"Holm-adjusted: Test 1 p = {adj[0]:.4g} -> {'PASS' if t1_pass else 'FAIL'}   "
      f"Test 2 p = {adj[1]:.4g} -> {'PASS' if t2_pass else 'FAIL'}")
print("-" * 72)

# ================= SECONDARY =================
print("\n" + "=" * 72)
print("SECONDARY (pre-declared, uncorrected)")
print("=" * 72)

sec = {}

# 8-gene panel
sec["panel8"] = {"auc_subtype": float(roc_auc_score(y_acc1, s8)),
                 "spearman_ttd": float(stats.spearmanr(s8[dead], td)[0]),
                 "spearman_p": float(stats.spearmanr(s8[dead], td)[1])}
print(f"8-gene panel        AUC={sec['panel8']['auc_subtype']:.3f}  "
      f"rho={sec['panel8']['spearman_ttd']:+.3f} (p={sec['panel8']['spearman_p']:.4g})")

# Ferrarotto MYC-TP63 (IN-SAMPLE comparator)
myc_tp63 = (Z["MYC"] - Z["TP63"]).values
sec["myc_tp63_insample"] = {"auc_subtype": float(roc_auc_score(y_acc1, myc_tp63)),
                            "spearman_ttd": float(stats.spearmanr(myc_tp63[dead], td)[0]),
                            "note": "IN-SAMPLE: derived on this cohort"}
print(f"MYC-TP63 [in-sample] AUC={sec['myc_tp63_insample']['auc_subtype']:.3f}  "
      f"rho={sec['myc_tp63_insample']['spearman_ttd']:+.3f}")

# Brayer 14 with published coefficients (IN-SAMPLE comparator)
bc = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_brayer/brayer_table5_coefs.csv")
b14 = bc[bc.in_14gene]
b14 = b14[b14.gene.isin(Z.columns)]
s_br = (Z[b14.gene.tolist()].values * b14.coefficient.values).sum(1)
auc_br = roc_auc_score(y_acc1, s_br)
if auc_br < 0.5:                      # orientation is not a free parameter for OUR panel,
    s_br, auc_br = -s_br, 1 - auc_br  # but their published score has no intercept
sec["brayer14_insample"] = {"auc_subtype": float(auc_br), "n_genes": int(len(b14)),
                            "spearman_ttd": float(stats.spearmanr(s_br[dead], td)[0]),
                            "note": "IN-SAMPLE: trained on these ACC-I/ACC-II labels"}
print(f"Brayer14 [in-sample] AUC={auc_br:.3f}  "
      f"rho={sec['brayer14_insample']['spearman_ttd']:+.3f}  ({len(b14)} genes)")

# convergence with their axis
r_myc = stats.spearmanr(s4, myc_tp63)
sec["ours_vs_myctp63"] = {"rho": float(r_myc[0]), "p": float(r_myc[1])}
print(f"\nours vs MYC-TP63    rho={r_myc[0]:+.3f} (p={r_myc[1]:.4g})")

# cross-tab of our median split vs their subtype
hi = (s4 > np.median(s4)).astype(int)
ct = pd.crosstab(pd.Series(hi, name="ours_high"), pd.Series(y_acc1, name="ACC_I"))
sec["crosstab_medsplit"] = ct.values.tolist()
sec["fisher_p"] = float(stats.fisher_exact(ct.values)[1])
print(f"median split vs ACC-I: Fisher p = {sec['fisher_p']:.4g}\n{ct}")

# adjusted for stage + solid histology
try:
    import statsmodels.api as sm
    solid = clin["Subtype"].astype(str).str.upper().eq("S").astype(float).values
    stg = clin["Stage"].astype(str).str.contains("IV|III").astype(float).values
    X = sm.add_constant(np.column_stack([s4, stg, solid]))
    m = sm.Logit(y_acc1, X).fit(disp=0)
    sec["adjusted_logit"] = {"score_coef": float(m.params[1]), "score_p": float(m.pvalues[1]),
                             "stage_p": float(m.pvalues[2]), "solid_p": float(m.pvalues[3])}
    print(f"adjusted logistic: score beta={m.params[1]:+.3f} p={m.pvalues[1]:.4g}  "
          f"(stage p={m.pvalues[2]:.3g}, solid p={m.pvalues[3]:.3g})")
except Exception as e:                                   # statsmodels optional
    sec["adjusted_logit"] = {"error": str(e)}
    print(f"adjusted logistic skipped: {e}")

# per-gene
sec["per_gene_auc"] = {g: float(roc_auc_score(y_acc1, Z[g].values)) for g in PANEL4}
print("per-gene AUC: " + "  ".join(f"{g}={v:.3f}" for g, v in sec["per_gene_auc"].items()))

res["secondary"] = sec

# ================= FIGURES =================
C1, C2 = "#B4436C", "#3C6E8F"

fig, ax = plt.subplots(1, 2, figsize=(9, 4))
for v, lab, c in [(s4[y_acc1 == 1], f"ACC-I (n={int(y_acc1.sum())})", C1),
                  (s4[y_acc1 == 0], f"ACC-II (n={int((1-y_acc1).sum())})", C2)]:
    ax[0].scatter(np.full(len(v), 0 if lab.startswith("ACC-I") else 1)
                  + RNG.normal(0, .06, len(v)), v, s=26, alpha=.8, color=c, label=lab,
                  edgecolor="none")
ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["ACC-I", "ACC-II"])
ax[0].set_ylabel("locked 4-gene score"); ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title(f"Test 1  AUC = {auc1:.3f} (p = {mw1.pvalue:.3g})", fontsize=10)

fpr = np.linspace(0, 1, 100)
o = np.argsort(-s4)
tp = np.cumsum(y_acc1[o]) / y_acc1.sum()
fp = np.cumsum(1 - y_acc1[o]) / (1 - y_acc1).sum()
ax[1].plot([0, 1], [0, 1], ls=":", c="0.6", lw=1)
ax[1].plot(fp, tp, c=C1, lw=2, label=f"4-gene (ours, external) {auc1:.3f}")
o8 = np.argsort(-myc_tp63)
ax[1].plot(np.cumsum(1 - y_acc1[o8]) / (1 - y_acc1).sum(), np.cumsum(y_acc1[o8]) / y_acc1.sum(),
           c=C2, lw=1.5, ls="--",
           label=f"MYC-TP63 (in-sample) {sec['myc_tp63_insample']['auc_subtype']:.3f}")
ax[1].set_xlabel("false positive rate"); ax[1].set_ylabel("true positive rate")
ax[1].legend(frameon=False, fontsize=8, loc="lower right")
ax[1].set_title("ROC, ACC-I vs ACC-II", fontsize=10)
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
save(fig, "CCR2020_test1_subtype")

fig, ax = plt.subplots(figsize=(4.6, 4))
ax.scatter(sd, td / 365.25, s=30, color=C1, alpha=.85, edgecolor="none")
b, a_ = np.polyfit(sd, td / 365.25, 1)
xs = np.linspace(sd.min(), sd.max(), 50)
ax.plot(xs, b * xs + a_, c="0.4", lw=1.2, ls="--")
ax.set_xlabel("locked 4-gene score"); ax.set_ylabel("time to death (years)")
ax.set_title(f"Test 2  Spearman rho = {rho:+.3f} (p = {p_rho:.3g}), n = {len(td)}", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
save(fig, "CCR2020_test2_timetodeath")

fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.hist(null_auc, bins=40, color="0.8", edgecolor="none")
ax.axvline(auc1, color=C1, lw=2)
ax.set_xlabel("AUC, random 4-gene panels"); ax.set_ylabel("count")
ax.set_title(f"Test 1 null: p = {p_rand:.4f}", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
save(fig, "CCR2020_null_random_panels")

# ---------------- write ----------------
pd.DataFrame({"sample": samples, "score4": s4, "score8": s8,
              "acc_subtype": subtype.values, "status": status,
              "ttdeath_days": ttd.values, "myc_tp63": myc_tp63}).to_csv(
    f"{OUT}/ccr2020_scores.csv", index=False)
json.dump(res, open(f"{OUT}/CCR2020_results.json", "w"), indent=1)
print(f"\nwrote {OUT}/CCR2020_results.json")

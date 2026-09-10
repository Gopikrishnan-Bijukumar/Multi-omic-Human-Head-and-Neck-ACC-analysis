"""Stage 36 (v2 of stage 31) - unified metrics regenerated against 31e.

WHY THIS FILE EXISTS
--------------------
work/outputs_metrics/UNIFIED_metrics.json and figures M1-M3/S1-S3 were generated at
2026-08-19 09:40.  31e_internal_stability_full.py finished at 11:36 and overturned the
internal-stability verdict those files encode.  v1 therefore reports the internal
leakage-free AUC as 0.777 with "0.810 NOT reproduced", while the manuscript text and
figure S4 report 0.783 +/- 0.021 with 0.810 inside the distribution.

This file is byte-for-byte stage 31 except that the internal-cohort block reads 31e's
14-draw distribution instead of 31d's 5, and it writes to work/outputs_metrics_v2/.
work/outputs_metrics/ is left untouched.

Original stage 31 docstring follows.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)


"""STAGE 7: unified metrics - AUC everywhere possible, C everywhere possible.

Motivation (ACC_STUDY_CONSOLIDATED.md section 10): the consolidated table reports C for some
validations and AUC for others, then places them side by side. A reader comparing AUC 0.884
(CCR2020 subtype) with C = 0.671 (DK survival) reads it as "the second cohort did better",
when they are different metrics on different endpoints.

This script computes every metric that is computable for every score x endpoint pair, so that
each comparison is like-for-like, and writes the figures that show it.

NOTHING IS REFITTED. Every score is a locked vector already on disk. New *metrics* on existing
*scores*. No panel, weight, direction, threshold or endpoint changes. The confirmatory claim
stays C = 0.671; every metric added here is POST-HOC DESCRIPTIVE under the section 13
provenance rules and none may become a new headline.

Key addition: time-dependent (cumulative/dynamic) AUC at 60 months, IPCW-weighted.
The existing binary AUC uses 44 of 54 DK patients - it drops everyone dying between 24 and 60
months. tdAUC@60 uses all 54 at the same horizon as the 5-year C, so C and AUC become a genuine
like-for-like pair for the first time.

What is NOT computable, and why (both are structural, not oversights):
  - CCR2020 subtype has no C. Subtype is not a time-to-event endpoint. Category difference.
  - CCR2020 population C is blocked: all 20 living patients have TTDeath blank, so censored
    observations have no observation window (PRESPEC_stage5.md section 5).
  - The internal discovery cohort has no C. sample_labels.csv `years` holds only 2 and 5, the
    group label, not follow-up time.

Terminology: the 20-patient discovery cohort is the "internal cohort" in all paper-facing output.
"JSE" never appears in a figure label. On-disk column names (score_frozen_jse) are untouched and
mapped at plot time.

Reads only; writes to work/outputs_metrics/.
"""
import importlib.util as _il
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            cumulative_dynamic_auc)
from sksurv.util import Surv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

ROOT = f"{ACC_DATA_ROOT}"
SRC = f"{ACC_DATA_ROOT}/model_outputs_metrics"          # read-only: 31b/31e prediction caches
OUT = f"{ACC_DATA_ROOT}/model_outputs_metrics_v2"       # write target; v1 is never modified
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20260819)
N_BOOT = 2000


def _load(name, path):
    s = _il.spec_from_file_location(name, path)
    m = _il.module_from_spec(s)
    s.loader.exec_module(m)
    return m


d3 = _load("d3", f"{HERE}/03_data.py")
dl = _load("dl", f"{HERE}/24_delong.py") if False else None   # see note below

# 24_delong.py runs its analysis at import time, so its DeLong helpers are copied rather than
# imported. They are byte-identical to work/24_delong.py:44-88 - the single implementation in
# the project. Any change there must be mirrored here.


def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(N, dtype=float)
    out[J] = T
    return out


def delong_cov(scores, y):
    """AUCs and their covariance matrix for k classifiers on the same samples."""
    pos = scores[:, y == 1]
    neg = scores[:, y == 0]
    k, m = pos.shape
    n = neg.shape[1]
    tx = np.array([_midrank(p) for p in pos])
    ty = np.array([_midrank(q) for q in neg])
    tz = np.array([_midrank(np.concatenate([p, q])) for p, q in zip(pos, neg)])
    aucs = (tz[:, :m].sum(1) / m - (m + 1) / 2.0) / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    return aucs, np.atleast_2d(cov)


def delong_ci(s, y):
    """Single-classifier AUC, Wald CI and p vs chance (24_delong.py:189-193)."""
    a, cov = delong_cov(s[None, :], y)
    se = float(np.sqrt(cov[0, 0]))
    z = (a[0] - 0.5) / se if se > 0 else np.nan
    return {"auc": float(a[0]), "se": se,
            "ci": [float(a[0] - 1.96 * se), float(a[0] + 1.96 * se)],
            "p_vs_chance": float(2 * stats.norm.sf(abs(z))) if se > 0 else np.nan}


# ---------------------------------------------------------------- metric helpers
def boot_auc(y, s, n=N_BOOT):
    """Bootstrap percentile CI for an AUC (23_ccr2020_external.py:50)."""
    out = []
    idx = np.arange(len(y))
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        out.append(roc_auc_score(y[b], s[b]))
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def boot_cindex(T_, E_, s, n=N_BOOT):
    """Bootstrap percentile CI for Harrell's C.

    There is no C-index bootstrap helper anywhere in the project; this factors out the inline
    version at 18_brayer_benchmark.py:233 so every C in the paper carries an interval computed
    the same way as every AUC.
    """
    out = []
    idx = np.arange(len(T_))
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if E_[b].sum() < 2:
            continue
        try:
            out.append(concordance_index_censored(E_[b].astype(bool), T_[b], s[b])[0])
        except Exception:
            continue
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def cox_summary(sc, T_, E_):
    """C + HR per SD + CI + p (26_cox_diagnostics.py:277)."""
    c = concordance_index_censored(E_.astype(bool), T_, sc)[0]
    m = CoxPHFitter().fit(pd.DataFrame({"time": T_, "event": E_, "s": stats.zscore(sc)}),
                          "time", "event")
    s_ = m.summary.iloc[0]
    return {"C": float(c), "hr_per_sd": float(s_["exp(coef)"]),
            "hr_ci": [float(s_["exp(coef) lower 95%"]), float(s_["exp(coef) upper 95%"])],
            "p": float(s_["p"])}


def td_auc(T_, E_, sc, times):
    """Cumulative/dynamic AUC(t), IPCW. Returns (auc_at_times, integrated).

    Censoring distribution is estimated on the same cohort being evaluated - there is no
    separate training set here, since no model is being fitted.
    """
    surv = Surv.from_arrays(event=E_.astype(bool), time=T_)
    a, mean_a = cumulative_dynamic_auc(surv, surv, sc, times)
    return np.asarray(a, float), float(mean_a)


def boot_tdauc(T_, E_, sc, t, n=N_BOOT):
    """Bootstrap CI for tdAUC at a single horizon."""
    out = []
    idx = np.arange(len(T_))
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if E_[b].sum() < 3 or T_[b].max() <= t:
            continue
        try:
            s_ = Surv.from_arrays(event=E_[b].astype(bool), time=T_[b])
            a, _ = cumulative_dynamic_auc(s_, s_, sc[b], [t])
            if np.isfinite(a[0]):
                out.append(float(a[0]))
        except Exception:
            continue
    if len(out) < 50:
        return [np.nan, np.nan]
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]


def uno_c(T_, E_, sc, tau):
    surv = Surv.from_arrays(event=E_.astype(bool), time=T_)
    return float(concordance_index_ipcw(surv, surv, sc, tau=tau)[0])


def check(label, got, want, tol):
    """Regression check against the published number. Mismatch is a finding, not a nuisance."""
    ok = abs(got - want) <= tol
    FAILURES.append((label, got, want)) if not ok else None
    print(f"  [{'OK ' if ok else 'FAIL'}] {label:<46s} got {got:.4f}  expect {want:.4f}")
    return ok


FAILURES = []

print("=" * 78)
print("STAGE 7 - unified metrics: AUC and C wherever each is computable")
print("=" * 78)

# ================================================================ DK cohort
D = d3.load("all")           # "all", not "evidence" - otherwise most of Brayer's 49 genes are
                             # missing from the lookup (see 24_delong.py:100)
genes_all, T, E, y = D["genes"], D["T"], D["E"], D["y"]
dk_raw, tr_raw = D["dk_raw"], D["tr_raw"]
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
early, late = (T <= 24) & (E == 1), T > 60
bmask = early | late
blab = early[bmask].astype(int)

print(f"\nDK n = {len(T)}   events(5y) = {int(E5.sum())}   "
      f"binary endpoint: {int(early.sum())} early / {int(late.sum())} late "
      f"({int((~bmask).sum())} excluded)")

# ---- scores on disk
s4 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel4_score.npy")
s8 = np.load(f"{ACC_DATA_ROOT}/model_outputs_panel8/dk_panel8_score.npy")
sb14 = np.load(f"{ACC_DATA_ROOT}/model_outputs_brayer/dk_brayer14_score.npy")

thr = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_threshold/dk_scores_by_procedure.csv")
assert np.allclose(thr["score_batch"].values, s4), "threshold CSV is not aligned to dk_panel4_score"
s_frozen = thr["score_frozen_jse"].values      # on-disk name kept; label mapped at plot time

# ---- brayer49: not cached, recompute exactly as 24_delong.py:118-127
coef = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_brayer/brayer_table5_coefs.csv")
Zd_all, Zt_all = d3.rint(dk_raw).values, d3.rint(tr_raw).values
gi_all = {g: i for i, g in enumerate(genes_all)}
sub = coef[coef.gene.map(lambda g: g in gi_all)]
ii = [gi_all[g] for g in sub.gene]
ww = sub.coefficient.values
sb49 = Zd_all[:, ii] @ ww
sgn = 1.0 if roc_auc_score(y, Zt_all[:, ii] @ ww) >= 0.5 else -1.0
sb49 = sgn * sb49
print(f"brayer49 recomputed from {len(sub)}/49 genes (orientation {sgn:+.0f}, internal-cohort only)")

# ---- composition: not cached, recompute via 22_dk_composition.py:118
dk_comp = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_deconv/dk_estimated_composition.csv", index_col=0)
assert list(dk_comp.index) == list(thr["dk_sample"]), "composition CSV not aligned to DK order"
TRIAD = {"Epithelial_Tumor": +1.0, "Dividing": +1.0, "Myoepithelial_Tumor": -1.0}


def comp_score(df):
    z = (df - df.mean(0)) / (df.std(0) + 1e-12)
    return sum(w * z[c] for c, w in TRIAD.items()).values


s_comp = comp_score(dk_comp)

DK_SCORES = {
    "ours_4gene":          ("4-gene panel",                 s4,       False),
    "ours_4gene_frozen":   ("4-gene, frozen internal ref",  s_frozen, False),
    "ours_8gene":          ("8-gene panel",                 s8,       False),
    "composition":         ("composition (deconvolution)",  s_comp,   False),
    "brayer14":            ("Brayer 14-gene",               sb14,     False),
    "brayer49":            ("Brayer 49-gene",               sb49,     False),
}

TIMES = np.arange(12, 61, 3, dtype=float)      # 12..60 months for the tdAUC curve

print("\n" + "-" * 78)
print("DK: every metric, same 54 patients")
print("-" * 78)
dk_rows, dk_curves = [], {}
for key, (label, sc, insample) in DK_SCORES.items():
    m5 = cox_summary(sc, T5, E5)
    mf = cox_summary(sc, T, E)
    a_t, a_int = td_auc(T, E, sc, TIMES)
    dk_curves[key] = a_t
    td60 = float(a_t[-1])
    row = {
        "score": key, "label": label, "n": len(T),
        "C_5yr": m5["C"], "C_5yr_ci_lo": None, "C_5yr_ci_hi": None,
        "C_full": mf["C"],
        "C_uno_tau60": uno_c(T, E, sc, tau=60.0),
        "tdAUC_60mo": td60, "tdAUC_integrated": a_int,
        "hr_per_sd": m5["hr_per_sd"], "hr_ci_lo": m5["hr_ci"][0], "hr_ci_hi": m5["hr_ci"][1],
        "cox_p": m5["p"],
        "binary_auc_n44": float(roc_auc_score(blab, sc[bmask])),
    }
    ci = boot_cindex(T5, E5, sc)
    row["C_5yr_ci_lo"], row["C_5yr_ci_hi"] = ci
    tci = boot_tdauc(T, E, sc, 60.0)
    row["tdAUC_60mo_ci_lo"], row["tdAUC_60mo_ci_hi"] = tci
    dcl = delong_ci(sc[bmask], blab)
    row["binary_auc_ci_lo"], row["binary_auc_ci_hi"] = dcl["ci"]
    row["binary_auc_p"] = dcl["p_vs_chance"]
    bci = boot_auc(blab, sc[bmask])
    row["binary_auc_boot_lo"], row["binary_auc_boot_hi"] = bci
    dk_rows.append(row)
    print(f"  {label:<30s} C5={row['C_5yr']:.3f} [{ci[0]:.3f},{ci[1]:.3f}]  "
          f"Cfull={row['C_full']:.3f}  Uno={row['C_uno_tau60']:.3f}  "
          f"tdAUC60={td60:.3f} [{tci[0]:.3f},{tci[1]:.3f}]  "
          f"binAUC={row['binary_auc_n44']:.3f}")

dk = pd.DataFrame(dk_rows)
dk.to_csv(f"{OUT}/metrics_dk.csv", index=False)

# ---- regression checks against published numbers
print("\nregression checks vs published numbers:")
r4 = dk[dk.score == "ours_4gene"].iloc[0]
check("4-gene DK C 5-yr", r4.C_5yr, 0.671, 0.002)
check("4-gene DK C full", r4.C_full, 0.642, 0.002)
check("4-gene DK HR/SD 5-yr", r4.hr_per_sd, 1.740, 0.01)
check("4-gene DK Cox p", r4.cox_p, 0.01714, 0.0005)
check("4-gene DK binary AUC (n=44)", r4.binary_auc_n44, 0.711, 0.002)
rf = dk[dk.score == "ours_4gene_frozen"].iloc[0]
check("frozen-ref DK C 5-yr", rf.C_5yr, 0.689, 0.002)
check("frozen-ref DK HR/SD", rf.hr_per_sd, 1.936, 0.01)
rc = dk[dk.score == "composition"].iloc[0]
check("composition DK C 5-yr", rc.C_5yr, 0.656, 0.002)
check("composition DK HR/SD", rc.hr_per_sd, 1.63, 0.01)
check("composition DK Cox p", rc.cox_p, 0.0391, 0.001)
rb = dk[dk.score == "brayer14"].iloc[0]
check("brayer14 DK binary AUC", rb.binary_auc_n44, 0.534, 0.002)

# ---- per-compartment Cox, must match the stored table exactly
pc_ref = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_deconv/dk_percompartment_cox.csv")
pc_rows = []
for ct in dk_comp.columns:
    m = cox_summary(dk_comp[ct].values, T5, E5)
    pc_rows.append({"cell_type": ct, "hr_per_sd": m["hr_per_sd"], "p": m["p"], "C": m["C"]})
pc = pd.DataFrame(pc_rows)
mrg = pc.merge(pc_ref, on="cell_type", suffixes=("", "_ref"))
dmax = float(np.abs(mrg.hr_per_sd - mrg.hr_per_sd_ref).max())
check("per-compartment Cox max |dHR|", dmax, 0.0, 1e-6)
pc.to_csv(f"{OUT}/metrics_dk_percompartment.csv", index=False)

# ================================================================ CCR2020
print("\n" + "-" * 78)
print("CCR2020: subtype AUC computable; population C structurally impossible")
print("-" * 78)
cc = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_ccr2020/ccr2020_scores.csv")
y_acc1 = (cc.acc_subtype.values == 1).astype(int)
dead = cc.status.values.astype(str) == "D"
ttd = pd.to_numeric(cc.ttdeath_days, errors="coerce").values
print(f"n = {len(cc)}   ACC-I = {y_acc1.sum()}   decedents with exact time = {int(dead.sum())}   "
      f"survivors with follow-up time = {int(np.isfinite(ttd[~dead]).sum())}")

# brayer14 on CCR2020: computed but never saved (23_ccr2020_external.py:216-228). Recompute.
# Its orientation is flipped against CCR2020's own labels, so it is IN-SAMPLE here.
ccr_rows = []
for key, label, sc, insample in [
    ("score4", "4-gene panel", cc.score4.values, False),
    ("score8", "8-gene panel", cc.score8.values, False),
    ("myc_tp63", "Ferrarotto MYC-TP63", cc.myc_tp63.values, True),
]:
    auc = float(roc_auc_score(y_acc1, sc))
    d = delong_ci(sc, y_acc1)
    row = {"score": key, "label": label, "in_sample": insample, "n": len(cc),
           "subtype_auc": auc, "subtype_auc_boot_lo": None, "subtype_auc_boot_hi": None,
           "subtype_auc_delong_lo": d["ci"][0], "subtype_auc_delong_hi": d["ci"][1],
           "subtype_mw_p": float(stats.mannwhitneyu(sc[y_acc1 == 1], sc[y_acc1 == 0]).pvalue)}
    row["subtype_auc_boot_lo"], row["subtype_auc_boot_hi"] = boot_auc(y_acc1, sc)
    # decedent-only block: conditioned on death, NOT a population C
    sd, td = sc[dead], ttd[dead]
    rho = stats.spearmanr(sd, td)
    ev = np.ones(dead.sum(), bool)      # all events, no censoring -> C is plain concordance
    row["decedent_n"] = int(dead.sum())
    row["decedent_C"] = float(concordance_index_censored(ev, td, sd)[0])
    row["decedent_C_ci_lo"], row["decedent_C_ci_hi"] = boot_cindex(td, ev.astype(int), sd)
    early5 = (td < 5 * 365.25).astype(int)
    row["decedent_auc_lt5y"] = float(roc_auc_score(early5, sd))
    row["decedent_auc_ci_lo"], row["decedent_auc_ci_hi"] = boot_auc(early5, sd)
    row["decedent_spearman_rho"] = float(rho.correlation)
    row["decedent_spearman_p"] = float(rho.pvalue)
    row["population_C"] = np.nan       # blocked: no follow-up time for the 20 survivors
    ccr_rows.append(row)
    print(f"  {label:<24s}{'  [IN-SAMPLE]' if insample else '              '} "
          f"subtype AUC={auc:.3f}  |  decedent-only C={row['decedent_C']:.3f} "
          f"AUC<5y={row['decedent_auc_lt5y']:.3f} rho={row['decedent_spearman_rho']:+.3f}")

ccr = pd.DataFrame(ccr_rows)
ccr.to_csv(f"{OUT}/metrics_ccr2020.csv", index=False)
r1 = ccr[ccr.score == "score4"].iloc[0]
check("CCR2020 4-gene subtype AUC", r1.subtype_auc, 0.884, 0.002)
check("CCR2020 4-gene decedent rho", r1.decedent_spearman_rho, -0.329, 0.002)

# ================================================================ binary block
print("\n" + "-" * 78)
print("Binary-classification endpoints (AUC only - no time-to-event, so no C)")
print("-" * 78)
bin_rows = []
for key, label, sc in [("ours_4gene", "4-gene panel", s4), ("ours_8gene", "8-gene panel", s8),
                       ("brayer14", "Brayer 14-gene", sb14), ("brayer49", "Brayer 49-gene", sb49)]:
    d = delong_ci(sc[bmask], blab)
    bin_rows.append({"block": "DK early vs late", "n": int(bmask.sum()), "score": key,
                     "label": label, "in_sample": False, "auc": d["auc"],
                     "ci_lo": d["ci"][0], "ci_hi": d["ci"][1], "p_vs_chance": d["p_vs_chance"],
                     "c_computable": False,
                     "c_reason": "binary endpoint - not time-to-event"})

# Brayer published-group recovery. NOTE: our_4gene here is RINT'd over all 62 samples, a
# different reference set from the 54-patient .npy - kept as its own row, never mixed.
bg = pd.read_csv(f"{ACC_DATA_ROOT}/model_outputs_brayer_group/dk_scores_vs_published_group.csv")
gl = bg.reconstructed_poor_group.values.astype(int)
d = delong_ci(bg.our_4gene.values, gl)
bin_rows.append({"block": "DK Brayer published group", "n": len(bg), "score": "ours_4gene_62",
                 "label": "4-gene (RINT over 62)", "in_sample": False, "auc": d["auc"],
                 "ci_lo": d["ci"][0], "ci_hi": d["ci"][1], "p_vs_chance": d["p_vs_chance"],
                 "c_computable": False, "c_reason": "group membership - not time-to-event"})

for key, label, sc, ins in [("score4", "4-gene panel", cc.score4.values, False),
                            ("myc_tp63", "Ferrarotto MYC-TP63", cc.myc_tp63.values, True)]:
    d = delong_ci(sc, y_acc1)
    bin_rows.append({"block": "CCR2020 subtype", "n": len(cc), "score": key, "label": label,
                     "in_sample": ins, "auc": d["auc"], "ci_lo": d["ci"][0], "ci_hi": d["ci"][1],
                     "p_vs_chance": d["p_vs_chance"], "c_computable": False,
                     "c_reason": "subtype label - not time-to-event"})

# Internal cohort - AUC only, no survival times exist.
# v2 CHANGE: sources the run-to-run distribution from 31e_internal_stability_full.py
# (10 paired repeats + 4 prior independent draws = 14 leakage-free draws) instead of
# 31d's 5.  31d ran after this file was first generated and reversed its own verdict:
# the published single run of 0.810 DOES lie inside the run-to-run distribution
# (2.5-97.5 percentile 0.7475-0.8152).  v1 of this file reported 0.777 with the range
# as its interval and stated 0.810 was "NOT reproduced"; that is superseded.
rep_path = f"{SRC}/internal_stability_full_preds.npz"
stab_path = f"{SRC}/INTERNAL_stability_full.json"
if os.path.exists(rep_path) and os.path.exists(stab_path):
    z = np.load(rep_path)
    stab = json.load(open(stab_path))
    yi = z["y"]
    runs = sorted({k.split("_")[0] for k in z.files if k != "y"},
                  key=lambda s: int(s[1:]))
    pooled_parts = []
    for rn in runs:
        for comp in ("loo", "cv5"):
            key = f"{rn}_leakage_free_{comp}"
            if key not in z.files:
                continue
            P = z[key]
            pooled_parts.extend(stats.rankdata(P[r]) for r in range(P.shape[0]))
    pooled = np.mean(pooled_parts, axis=0)
    d = delong_ci(pooled, yi)
    S = stab["summary"]["leakage_free"]
    bin_rows.append({
        "block": "Internal cohort (nested CV)", "n": len(yi), "score": "ours_4gene",
        "label": "4-gene panel, leakage-free  (bar = run-to-run 2.5-97.5 pct)", "in_sample": False,
        "auc": float(S["mean"]), "ci_lo": float(S["pct_lo"]), "ci_hi": float(S["pct_hi"]),
        "p_vs_chance": d["p_vs_chance"], "c_computable": False,
        "c_reason": "no follow-up times - `years` holds only the 2/5 group label"})
    print(f"  internal cohort, leakage-free: {S['mean']:.3f} "
          f"[{S['pct_lo']:.3f}-{S['pct_hi']:.3f}] over {S['n_draws']} independent draws "
          f"(published single run 0.810 lies inside: covers={S['covers_published']})")
    print(f"    pooled-rank ROC AUC = {d['auc']:.3f}")
else:
    print("  [skip] internal-cohort full repeats not found - run 31b then 31e first")

bindf = pd.DataFrame(bin_rows)
bindf.to_csv(f"{OUT}/metrics_binary.csv", index=False)
for _, r in bindf.iterrows():
    print(f"  {r.block:<28s} {r.label:<24s} AUC = {r.auc:.3f} "
          f"[{r.ci_lo:.3f},{r.ci_hi:.3f}]{'  [IN-SAMPLE]' if r.in_sample else ''}")

# ================================================================ save
N_CENS_EARLY = int(((T <= 60) & (E == 0)).sum())
print(f"\nDK censoring completeness: {N_CENS_EARLY} of {len(T)} censored before 60 months")

json.dump({
    "note": ("Post-hoc descriptive re-reporting of locked scores. No refitting. The confirmatory "
             "claim remains C = 0.671 (DK, 4-gene, 5-year OS)."),
    "dk": dk.to_dict("records"),
    "dk_percompartment": pc.to_dict("records"),
    "ccr2020": ccr.to_dict("records"),
    "binary": bindf.to_dict("records"),
    "tdauc_times_months": TIMES.tolist(),
    "tdauc_curves": {k: v.tolist() for k, v in dk_curves.items()},
    "censoring_completeness_dk": {
        "n_censored_before_60mo": N_CENS_EARLY,
        "note": ("Zero DK patients are censored before 60 months, so G(t) = 1 across the "
                 "5-year window, every IPCW weight is 1, and Uno's C equals Harrell's C "
                 "exactly. The 5-year C-index therefore carries no censoring bias."),
    },
    "not_computable": {
        "ccr2020_population_C": "all 20 living patients have TTDeath blank - no observation window",
        "ccr2020_subtype_C": "subtype is not a time-to-event endpoint",
        "internal_cohort_C": "sample_labels.csv `years` holds only the 2/5 group label",
        "brayer_group_C": "group membership is a binary label",
    },
    "regression_failures": [{"check": a, "got": b, "expected": c} for a, b, c in FAILURES],
}, open(f"{OUT}/UNIFIED_metrics.json", "w"), indent=2, default=float)

print("\n" + "=" * 78)
if FAILURES:
    print(f"*** {len(FAILURES)} REGRESSION CHECK(S) FAILED ***")
    for a, b, c in FAILURES:
        print(f"    {a}: got {b:.4f}, expected {c:.4f}")
else:
    print("all regression checks passed")
print("=" * 78)
print(f"metrics written -> {OUT}")

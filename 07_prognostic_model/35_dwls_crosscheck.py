"""Stage 35 - cross-check the local DWLS reimplementation against the published package.

The review's deconvolution section says the implementation "is a local reimplementation of
dampened weighted least squares rather than the published package and should be cross-checked
before submission". At the time of the review no R deconvolution package was installed here.
DWLS 0.1.0 is current on CRAN (not archived, as work/README and the consolidated document
assumed), so this stage uses the published solveDampenedWLS directly.

Three tests:
  1. Agreement - same signature, same bulk, published solver vs work/22_dk_composition.py::dwls.
  2. Consequence - does the composition survival result (C = 0.656, HR/SD = 1.63, p = 0.039)
     survive when the reference implementation supplies the proportions?
  3. Ground truth - leave-one-patient-out synthetic mixtures at KNOWN proportions, scoring both
     implementations against truth. This tests the solver directly rather than through the
     internal-cohort gate, whose absolute-composition criterion (G3) failed at 0.685 vs 0.70.

Writes to work/outputs_dwls/. work/outputs_deconv/ is read-only here.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import os, sys, json, subprocess
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import nnls
from sksurv.metrics import concordance_index_censored
from lifelines import CoxPHFitter
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
SRC22 = f"{HERE}/22_dk_composition.py"

# Reuse stage 22's own signature builder and solver verbatim (the 25b/25c pattern), so the
# comparison is against the exact code that produced the published composition numbers.
_src = open(SRC22).read().split("# ---- survival + expression score (read-only) ----")[0]
exec(compile(_src, SRC22, "exec"))          # defines COARSE, pb, build_signature, dwls, S, shared, dk_log, tr_log

OUT = f"{ACC_DATA_ROOT}/model_outputs_dwls"               # override stage 22's OUT before anything is written
os.makedirs(OUT, exist_ok=True)
DECONV = f"{ACC_DATA_ROOT}/model_outputs_deconv"
rng = np.random.default_rng(35)
TRIAD = {"Epithelial_Tumor": +1.0, "Dividing": +1.0, "Myoepithelial_Tumor": -1.0}
RSCRIPT = f"{HERE}/35r_dwls_reference.R"

import importlib.util
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

D = d3.load("all")
T, E = D["T"], D["E"]
dk_raw = D["dk_raw"]
dk_ids = list(dk_raw.index)
T5, E5 = np.minimum(T, 60.0), np.where(T > 60, 0, E)
dk_lin = (2.0 ** dk_log.loc[dk_ids, shared]) - 1.0; dk_lin[dk_lin < 0] = 0
print(f"\nDK n={len(dk_ids)} events5y={E5.sum()}  signature {S.shape[0]}x{S.shape[1]}", flush=True)

res = {"_note": "Stage 35 DWLS cross-check. outputs_deconv/ is read-only here.",
       "dwls_package": "CRAN DWLS 0.1.0 (solveDampenedWLS)",
       "signature": {"n_markers": int(S.shape[0]), "compartments": list(S.columns)}}


def run_reference(sig, bulk, tag):
    """bulk: DataFrame genes x samples."""
    sf, bf, of = f"{OUT}/_sig_{tag}.csv", f"{OUT}/_bulk_{tag}.csv", f"{OUT}/_ref_{tag}.csv"
    sig.to_csv(sf); bulk.to_csv(bf)
    p = subprocess.run(["Rscript", RSCRIPT, sf, bf, of], capture_output=True, text=True)
    print(p.stdout.strip()[-400:] or p.stderr.strip()[-400:], flush=True)
    if not os.path.exists(of):
        raise RuntimeError(f"reference DWLS failed for {tag}:\n{p.stderr[-2000:]}")
    out = pd.read_csv(of, index_col=0)
    for f in (sf, bf):
        os.remove(f)
    return out


def comp_score(df):
    z = (df - df.mean(0)) / (df.std(0) + 1e-12)
    return sum(w * z[c] for c, w in TRIAD.items()).values


def surv_stats(sc, mask=None):
    t5, e5 = (T5, E5) if mask is None else (T5[mask], E5[mask])
    c = float(concordance_index_censored(e5.astype(bool), t5, sc)[0])
    cph = CoxPHFitter().fit(pd.DataFrame({"time": t5, "event": e5, "s": stats.zscore(sc)}),
                            "time", "event")
    return {"C_5yr": c, "hr_per_sd": float(np.exp(cph.params_.iloc[0])),
            "cox_p": float(cph.summary["p"].iloc[0])}


# =============================================================== 1. agreement on DK
print("\n[1] agreement: published solveDampenedWLS vs local dwls(), same signature/bulk", flush=True)
local = pd.DataFrame({p: dwls(dk_lin.loc[p, S.index], S) for p in dk_ids}).T
stored = pd.read_csv(f"{DECONV}/dk_estimated_composition.csv", index_col=0)
repro = float(np.max(np.abs(local[stored.columns].values - stored.values)))
print(f"  local solver reproduces stored composition to {repro:.2e}", flush=True)

ref = run_reference(S, dk_lin[S.index].T, "dk")
ref = ref[local.columns]
ref.to_csv(f"{OUT}/dk_composition_reference_dwls.csv")
local.to_csv(f"{OUT}/dk_composition_local_dwls.csv")

# The published solver's internal quadprog can fail on individual samples
# ("constraints are inconsistent"). Those rows come back all-NaN. Record them and
# restrict every comparison below to the samples both implementations solved.
failed = [i for i in ref.index if ref.loc[i].isna().all()]
ok = [i for i in local.index if i not in failed]
ok_mask = np.array([i in set(ok) for i in dk_ids])
print(f"  reference solver failed on {len(failed)}/{len(local)} samples: {failed}", flush=True)

per_ct = []
for c in local.columns:
    a, b = local.loc[ok, c].values, ref.loc[ok, c].values
    per_ct.append({"compartment": c, "mean_local": float(a.mean()), "mean_reference": float(b.mean()),
                   "pearson": float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan,
                   "spearman": float(stats.spearmanr(a, b).statistic),
                   "max_abs_diff": float(np.max(np.abs(a - b))),
                   "mean_abs_diff": float(np.mean(np.abs(a - b)))})
per_ct = pd.DataFrame(per_ct).sort_values("mean_local", ascending=False)
per_ct.to_csv(f"{OUT}/agreement_per_compartment.csv", index=False)
print(per_ct.to_string(index=False), flush=True)
res["agreement_dk"] = {
    "local_reproduces_stored_max_abs_diff": repro,
    "reference_solver_failures": {"n": len(failed), "samples": failed,
                                  "reason": "quadprog inside DWLS::solveDampenedWLS reported "
                                            "'constraints are inconsistent, no solution!'",
                                  "n_compared": len(ok)},
    "per_compartment": per_ct.to_dict("records"),
    "median_spearman": float(per_ct.spearman.median()),
    "max_abs_diff_overall": float(per_ct.max_abs_diff.max()),
    "triad_spearman": {c: float(per_ct.set_index("compartment").loc[c, "spearman"]) for c in TRIAD}}

# =============================================================== 2. does the result survive?
print("\n[2] consequence for the published composition survival result", flush=True)
# composition scores are z-scored within cohort, so both are computed on the same
# solved subset to keep the comparison like-for-like
sc_local_all = comp_score(local)
sc_local = comp_score(local.loc[ok])
sc_ref = comp_score(ref.loc[ok])
st_local_all = surv_stats(sc_local_all, None)
st_local = surv_stats(sc_local, ok_mask)
st_ref = surv_stats(sc_ref, ok_mask)
published = json.load(open(f"{DECONV}/DK_COMPOSITION_results.json"))
res["composition_survival"] = {
    "n_compared": len(ok),
    "local_all_54": st_local_all,
    "local": st_local, "reference": st_ref,
    "score_spearman_local_vs_reference": float(stats.spearmanr(sc_local, sc_ref).statistic),
    "published_reference_values": {"C_5yr": 0.65590, "hr_per_sd": 1.63027, "cox_p": 0.03906}}
print(f"  local(54) C={st_local_all['C_5yr']:.4f} HR/SD={st_local_all['hr_per_sd']:.3f} p={st_local_all['cox_p']:.4f}", flush=True)
print(f"  local     C={st_local['C_5yr']:.4f} HR/SD={st_local['hr_per_sd']:.3f} p={st_local['cox_p']:.4f}", flush=True)
print(f"  reference C={st_ref['C_5yr']:.4f} HR/SD={st_ref['hr_per_sd']:.3f} p={st_ref['cox_p']:.4f}", flush=True)
print(f"  score rank agreement rho={res['composition_survival']['score_spearman_local_vs_reference']:.4f}", flush=True)

# =============================================================== 3. synthetic ground truth
print("\n[3] leave-one-patient-out synthetic mixtures at known proportions", flush=True)
pats = sorted(pb.index.get_level_values("patient").unique())
cts = sorted(pb.index.get_level_values("cell_type").unique())


def signature_excluding(pat, n_markers=60):
    sub = pb[pb.index.get_level_values("patient") != pat]
    prof = sub.groupby(level="cell_type").sum()
    prof = prof.div(prof.sum(axis=1), axis=0) * 1e6
    lv = np.log2(prof + 1.0)
    markers = []
    for ct in prof.index:
        spec = lv.loc[ct] - lv.drop(index=ct).max(axis=0)
        markers += list(spec[prof.loc[ct] > 20].sort_values(ascending=False).index[:n_markers])
    return prof.T.loc[sorted(set(markers))]


N_MIX = 4
truth_rows, mix_cols, mix_sig = [], {}, {}
for pat in pats:
    prof = pb.xs(pat, level="patient")
    prof = prof.reindex(cts).dropna(how="all")
    if prof.shape[0] < 3:
        continue
    prof = prof.div(prof.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0) * 1e6
    Sx = signature_excluding(pat)
    for m in range(N_MIX):
        w = rng.dirichlet(np.ones(prof.shape[0]) * 0.7)
        mix = (prof.T.values @ w)
        name = f"{pat}__m{m}"
        mix_cols[name] = pd.Series(mix, index=prof.columns)
        mix_sig[name] = Sx
        truth_rows.append(dict({"sample": name, "patient": pat},
                               **{ct: float(v) for ct, v in zip(prof.index, w)}))
truth = pd.DataFrame(truth_rows).set_index("sample").fillna(0.0)
print(f"  {len(truth)} synthetic mixtures from {truth.patient.nunique()} patients", flush=True)

est_local, est_ref = {}, {}
for pat in truth.patient.unique():
    names = truth.index[truth.patient == pat]
    Sx = mix_sig[names[0]]
    Bx = pd.DataFrame({n: mix_cols[n] for n in names}).reindex(Sx.index).fillna(0.0)
    for n in names:
        est_local[n] = dwls(Bx[n], Sx)
    rr = run_reference(Sx, Bx, f"syn_{pat}")
    for n in names:
        if n in rr.index and not rr.loc[n].isna().all():
            est_ref[n] = rr.loc[n]
EL = pd.DataFrame(est_local).T
ER = pd.DataFrame(est_ref).T
syn_failed = [i for i in EL.index if i not in ER.index]
if syn_failed:
    print(f"  reference solver failed on {len(syn_failed)} synthetic mixtures", flush=True)
keep_syn = [i for i in EL.index if i in ER.index]
EL, ER, truth = EL.loc[keep_syn], ER.loc[keep_syn], truth.loc[keep_syn]
cols = [c for c in truth.columns if c != "patient"]
EL = EL.reindex(columns=cols).fillna(0.0); ER = ER.reindex(columns=cols).fillna(0.0)
TT = truth[cols]
EL.to_csv(f"{OUT}/synthetic_estimate_local.csv"); ER.to_csv(f"{OUT}/synthetic_estimate_reference.csv")
TT.to_csv(f"{OUT}/synthetic_truth.csv")


def recovery(est):
    per_sample = [float(stats.spearmanr(TT.loc[i], est.loc[i]).statistic) for i in TT.index]
    per_ct = {c: float(stats.spearmanr(TT[c], est[c]).statistic) for c in cols if TT[c].std() > 0}
    rmse = float(np.sqrt(((TT.values - est.values) ** 2).mean()))
    return {"mean_per_sample_spearman": float(np.nanmean(per_sample)),
            "median_per_sample_spearman": float(np.nanmedian(per_sample)),
            "per_compartment_spearman": per_ct,
            "median_per_compartment_spearman": float(np.nanmedian(list(per_ct.values()))),
            "rmse": rmse,
            "mean_abs_error": float(np.abs(TT.values - est.values).mean()),
            "triad_bias": {c: float((est[c] - TT[c]).mean()) for c in TRIAD if c in cols}}


res["synthetic_recovery"] = {"n_mixtures": int(len(TT)),
                             "reference_solver_failures": len(syn_failed), "n_patients": int(truth.patient.nunique()),
                             "design": "LOO signature: the mixed patient never contributes to its own reference",
                             "local": recovery(EL), "reference": recovery(ER),
                             "local_vs_reference_mean_abs_diff": float(np.abs(EL.values - ER.values).mean())}
for k in ("local", "reference"):
    d = res["synthetic_recovery"][k]
    print(f"  {k:9s} per-sample rho {d['mean_per_sample_spearman']:.3f} | "
          f"per-compartment rho {d['median_per_compartment_spearman']:.3f} | RMSE {d['rmse']:.4f}", flush=True)

for f in os.listdir(OUT):
    if f.startswith("_ref_") or f.startswith("_sig_") or f.startswith("_bulk_"):
        os.remove(f"{OUT}/{f}")
json.dump(res, open(f"{OUT}/DWLS_CROSSCHECK.json", "w"), indent=2)
print(f"\nwrote {OUT}/DWLS_CROSSCHECK.json", flush=True)

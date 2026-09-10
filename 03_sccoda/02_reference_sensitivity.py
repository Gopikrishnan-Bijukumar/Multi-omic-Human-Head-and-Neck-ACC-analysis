import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import os, sys, json, time
os.environ["TF_CPP_MIN_LOG_LEVEL"]="3"
import anndata as ad, pandas as pd, numpy as np, tensorflow as tf
from sccoda.util import cell_composition_data as dat, comp_ana as mod

A = ad.read_h5ad(f"{ACC_DATA_ROOT}/outputs/h5ad_files/inprogress_3.h5ad", backed="r")
obs = A.obs
META = obs.groupby("batch", observed=True)[["clinical_outcome","model_outcome","site"]].first()

def batch_means_se(ind, nb=50):
    m = len(ind)//nb
    bm = np.array([ind[i*m:(i+1)*m].mean() for i in range(nb)])
    return bm.std(ddof=1)/np.sqrt(nb)

def run(tag, level, ref, seed, drop_batches=(), outcome="clinical_outcome", nres=100_000, nburn=10_000):
    o = obs[~obs["batch"].isin(drop_batches)]
    comp = pd.crosstab(o["batch"], o[level])
    comp = comp.loc[:, comp.sum(0) > 0]
    meta = META.reindex(comp.index)
    y = pd.Categorical(meta[outcome].astype(str), categories=["Good","Poor"], ordered=True)
    cdf = comp.copy(); cdf[outcome] = y
    d = dat.from_pandas(cdf, covariate_columns=[outcome])
    np.random.seed(seed); tf.random.set_seed(seed)
    m = mod.CompositionalAnalysis(d, formula=outcome, reference_cell_type=ref)
    refname = str(d.var_names[int(m.reference_cell_type)])
    t0=time.time(); r = m.sample_hmc(num_results=nres, num_burnin=nburn, verbose=False); el=time.time()-t0
    _, eff = r.summary_prepare(est_fdr=0.05)
    thr = r.model_specs["threshold_prob"]
    beta = np.asarray(r.posterior["beta"])[0]
    ses={}
    for k, ct in enumerate(d.var_names):
        ind = (np.abs(beta[:,0,k])>1e-3).astype(float)
        ses[str(ct)] = batch_means_se(ind)
    e = eff.reset_index()
    e["mcse_inclusion"]=[ses[str(c)] for c in e["Cell Type"]]
    e["credible"]=e["Final Parameter"].ne(0)
    e["tag"]=tag; e["level"]=level; e["reference"]=refname; e["seed"]=seed
    e["threshold"]=thr; e["n_good"]=int((y=="Good").sum()); e["n_poor"]=int((y=="Poor").sum())
    cred=e.loc[e.credible,"Cell Type"].tolist()
    print(f"\n===== {tag} | level={level} ref={refname} seed={seed} n={len(comp)} ({int((y=='Good').sum())}G/{int((y=='Poor').sum())}P) thr={thr:.3f} {el:.0f}s",flush=True)
    top=e.sort_values("Inclusion probability",ascending=False).head(6)
    for _,row in top.iterrows():
        print(f"   {row['Cell Type'][:45]:47s} incl={row['Inclusion probability']:.4f} (+-{row['mcse_inclusion']:.4f})  logFC={row['log2-fold change']:+.3f}  {'CREDIBLE' if row['credible'] else ''}",flush=True)
    print("   CREDIBLE SET:", cred, flush=True)
    return e

MAXILLA = META.index[META["site"].astype(str)=="Maxilla"].tolist()
print("Maxilla patients (all Poor):", MAXILLA, flush=True)
runs=[]
runs.append(run("A5_mural_s0","annot_5","automatic",20260729))
runs.append(run("A5_mural_s1","annot_5","automatic",1))
runs.append(run("A5_mural_s2","annot_5","automatic",2))
runs.append(run("A5_endo","annot_5","Endothelial cells",20260729))
runs.append(run("A5_m2","annot_5","Macrophage - M2",20260729))
runs.append(run("A5_fibro","annot_5","Fibroblast cells",20260729))
runs.append(run("A6_mural_s0","annot_6","automatic",20260729))
runs.append(run("A6_mural_s1","annot_6","automatic",1))
runs.append(run("A6_mural_s2","annot_6","automatic",2))
runs.append(run("A5_noMaxilla","annot_5","automatic",20260729,drop_batches=MAXILLA))
runs.append(run("A6_noMaxilla","annot_6","automatic",20260729,drop_batches=MAXILLA))
out=pd.concat(runs,ignore_index=True)
out.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"sensitivity_results.csv"),index=False)
print("\nDONE -> sensitivity_results.csv",flush=True)

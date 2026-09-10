"""Which biological axis actually carries prognosis in DK? Fixed, hand-specified gene
sets (no fitting at all) scored two ways, plus WGCNA module eigengenes projected into DK.

Nothing is fitted to DK survival; these are pre-specified sets scored once.
"""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

from pathlib import Path as _Path
# Directory holding this script, so sibling modules resolve wherever the repo lives.
HERE = str(_Path(__file__).resolve().parent)

import sys, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sksurv.metrics import concordance_index_censored
import warnings; warnings.filterwarnings("ignore")

R = f"{ACC_DATA_ROOT}"
spec = importlib.util.spec_from_file_location("d3", f"{HERE}/03_data.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)

D = d3.load("all")
dk_raw, tr_raw, T, E = D["dk_raw"], D["tr_raw"], D["T"], D["E"]
Zd = d3.rint(dk_raw)
print("DK", Zd.shape, "events", E.sum())

SETS = {
 "cell_cycle_G2M": "AURKA AURKB BUB1 BUB1B CCNB1 CCNB2 CDC20 CDCA8 CDK1 CENPA CENPE CENPF KIF11 KIF20A KIF23 MKI67 NDC80 NUSAP1 PLK1 PRC1 RRM2 TOP2A TPX2 TTK UBE2C ANLN PIMREG NCAPG",
 "cell_cycle_S": "MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 PCNA RRM1 TYMS GINS1 GINS2 CDC45 CDT1 EXO1 FEN1 RFC4 POLA1 CHAF1A CLSPN DSCC1",
 "E2F_targets": "E2F1 E2F2 MYBL2 CCNE1 CCNE2 CDC6 ORC1 ORC6 DHFR TK1 CDKN2A FOXM1 MCM10 BRCA1 BRCA2 RAD51 CHEK1",
 "EMT": "VIM FN1 CDH2 SNAI1 SNAI2 TWIST1 ZEB1 ZEB2 SPARC COL1A1 COL1A2 COL5A1 TIMP1 SERPINE1 TAGLN TPM1 THBS2 POSTN LOXL2 ITGB1",
 "hypoxia": "CA9 VEGFA SLC2A1 LDHA PGK1 ENO1 ALDOA ADM P4HA1 P4HA2 NDRG1 BNIP3 EGLN3 HK2 PDK1 ANKRD37",
 "myoepithelial": "TP63 KRT5 KRT14 KRT17 ACTA2 MYH11 CNN1 TAGLN OXTR MYLK CAV1 SPARCL1 DST LAMB3 ITGA6",
 "luminal_epithelial": "KIT KRT7 KRT8 KRT18 CD24 EPCAM MUC1 SOX10 ELF5 AR",
 "MYB_program": "MYB MYBL1 NFIB KIT EN1 SOX4 CCNB1 BCL2 MYC RUNX1 CDK6",
 "NOTCH_activation": "HES1 HES4 HEY1 HEY2 NOTCH1 NOTCH3 DTX1 NRARP JAG1 DLL1 DLL3 LFNG",
 "immune_cytotoxic": "CD8A CD3D CD3E GZMB PRF1 IFNG CXCL9 CXCL10 GZMK LCK ZAP70 CD2 IL2RB",
 "immune_myeloid": "CD68 CD163 CSF1R AIF1 ITGAM MRC1 TYROBP FCER1G LYZ C1QA C1QB",
 "stress_ISG": "ISG15 IFI6 IFI44L MX1 OAS1 OAS2 STAT1 IRF7 HERC5 IFIT1 IFIT3 RSAD2",
 "hedgehog_stem": "SOX2 POU5F1 NANOG PROM1 ALDH1A1 CD44 BMI1 NOTCH1",
 "neuronal_perineural": "NCAM1 L1CAM GFRA1 GFRA3 NTRK2 NGFR S100B PLP1 MPZ SOX10 GAP43",
 "translation_MYC": "MYC NPM1 NCL RPL3 RPS6 EIF4E EIF4G1 EIF4A1 PA2G4 SRM ODC1 TFAP4 MAGOH",
 "TERT_telomere": "TERT TERF1 TERF2 POT1 RTEL1 DKC1",
 "AXL_receptor": "AXL GAS6 MERTK TYRO3",
 "senescence_SASP": "CDKN1A CDKN2A IL6 IL1B CXCL8 MMP3 MMP9 SERPINE1 IGFBP3",
}


def score_rint(genes):
    g = [x for x in genes if x in Zd.columns]
    return Zd[g].mean(axis=1).values, len(g)


def score_singscore(genes):
    """Within-sample rank score: normalisation-free, the cross-platform standard."""
    ranks = dk_raw.rank(axis=1, pct=True)
    g = [x for x in genes if x in ranks.columns]
    return ranks[g].mean(axis=1).values, len(g)


def ev(sc, name, extra=""):
    if np.std(sc) < 1e-12: return None
    c = concordance_index_censored(E.astype(bool), T, sc)[0]
    df = pd.DataFrame({"time": T, "event": E, "score": stats.zscore(sc)})
    cph = CoxPHFitter().fit(df, "time", "event")
    hi = sc > np.median(sc)
    lr = logrank_test(T[hi], T[~hi], E[hi], E[~hi])
    return dict(sig=name, n_genes=extra, C=c, hr_per_sd=np.exp(cph.params_.iloc[0]),
                cox_p=cph.summary["p"].iloc[0], logrank_p=lr.p_value)


rows = []
for name, gs in SETS.items():
    genes = gs.split()
    for tag, fn in [("rint", score_rint), ("singscore", score_singscore)]:
        sc, n = fn(genes)
        r = ev(sc, f"{name}[{tag}]", n)
        if r: rows.append(r)

# ---- WGCNA modules projected into DK (module eigengene = PC1 of module genes) ----
wg = pd.read_csv(f"{ACC_DATA_ROOT}/bulk_rna_data/outputs/images/clinical_outcome/wgcna/06_gene_info_complete.csv")
gcol = [c for c in wg.columns if c.lower() in ("gene", "genes", "symbol", "gene_name", "x")]
print("wgcna cols:", wg.columns[:8].tolist(), "modules:", wg["Module"].value_counts().to_dict() if "Module" in wg else None)
gc = gcol[0] if gcol else wg.columns[0]
for mod, sub in wg.groupby("Module"):
    genes = [g for g in sub[gc].astype(str) if g in Zd.columns]
    if len(genes) < 15: continue
    X = Zd[genes].values
    u, s_, vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    pc1 = u[:, 0] * s_[0]
    # orient by module membership sign convention: correlate with mean expression
    if np.corrcoef(pc1, X.mean(1))[0, 1] < 0: pc1 = -pc1
    r = ev(pc1, f"WGCNA_{mod}_eigengene", len(genes))
    if r: rows.append(r)

res = pd.DataFrame(rows).sort_values("cox_p")
pd.set_option("display.width", 200)
print(res.to_string(index=False))
res.to_csv(f"{ACC_DATA_ROOT}/model_outputs/biology_axis_scan.csv", index=False)

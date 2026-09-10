"""Module design matrix: the middle level of the multi-level model.

Modules come from three DK-independent sources:
  1. curated biological programs (hallmark-style, hand-specified)
  2. WGCNA modules from the discovery bulk cohort
  3. single-cell pseudobulk DE gene sets (per enriched cell type, per direction)
  4. the evidence-consensus direction sets (Poor-up / Good-up)

None of them is derived from DK survival.
"""
import numpy as np, pandas as pd

CURATED = {
 "hypoxia": "CA9 VEGFA SLC2A1 LDHA PGK1 ENO1 ALDOA ADM P4HA1 P4HA2 NDRG1 BNIP3 EGLN3 HK2 PDK1 ANKRD37 PFKFB3 PFKFB4 ALDOC BNIP3L HILPDA VHL EGLN1 CAV1 ANGPTL4 SLC16A3 MIF DDIT4",
 "glycolysis": "HK1 HK2 GPI PFKL PFKP ALDOA TPI1 GAPDH PGK1 PGAM1 ENO1 ENO2 PKM LDHA SLC2A1 SLC16A1 SLC16A3",
 "cell_cycle_G2M": "AURKA AURKB BUB1 BUB1B CCNB1 CCNB2 CDC20 CDCA8 CDK1 CENPA CENPE CENPF KIF11 KIF20A KIF23 MKI67 NDC80 NUSAP1 PLK1 PRC1 RRM2 TOP2A TPX2 TTK UBE2C ANLN PIMREG NCAPG NCAPG2 KIF2C",
 "cell_cycle_S": "MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 PCNA RRM1 TYMS GINS1 GINS2 CDC45 CDT1 EXO1 FEN1 RFC4 POLA1 CHAF1A CLSPN DSCC1 PRIM1 POLE2",
 "E2F_targets": "E2F1 E2F2 MYBL2 CCNE1 CCNE2 CDC6 ORC1 ORC6 DHFR TK1 FOXM1 MCM10 BRCA1 BRCA2 RAD51 CHEK1 CDC25A TIPIN",
 "MYC_translation": "MYC NPM1 NCL RPL3 RPS6 EIF4E EIF4G1 EIF4A1 PA2G4 SRM ODC1 TFAP4 MAGOH NOP56 NOP58 FBL DKC1 GNL3 POLR1B RRP9 BYSL",
 "EMT": "VIM FN1 CDH2 SNAI1 SNAI2 TWIST1 ZEB1 ZEB2 SPARC COL1A1 COL1A2 COL5A1 TIMP1 SERPINE1 TAGLN TPM1 THBS2 POSTN LOXL2 ITGB1 LOX FBN1 COL6A3",
 "myoepithelial": "TP63 KRT5 KRT14 KRT17 ACTA2 MYH11 CNN1 TAGLN OXTR MYLK CAV1 SPARCL1 DST LAMB3 ITGA6 MYL9 CAV2",
 "luminal_epithelial": "KIT KRT7 KRT8 KRT18 CD24 EPCAM MUC1 SOX10 ELF5 AR CLDN3 CLDN4 CLDN7",
 "MYB_program": "MYB MYBL1 NFIB KIT EN1 SOX4 BCL2 RUNX1 CDK6 CCND1",
 "NOTCH": "HES1 HES4 HEY1 HEY2 NOTCH1 NOTCH3 DTX1 NRARP JAG1 DLL1 DLL3 LFNG MAML2",
 "immune_cytotoxic": "CD8A CD3D CD3E GZMB PRF1 IFNG CXCL9 CXCL10 GZMK LCK ZAP70 CD2 IL2RB CD247 GZMA",
 "immune_myeloid": "CD68 CD163 CSF1R AIF1 ITGAM MRC1 TYROBP FCER1G LYZ C1QA C1QB C1QC MSR1",
 "interferon": "ISG15 IFI6 IFI44L MX1 MX2 OAS1 OAS2 OAS3 STAT1 STAT2 IRF7 IFIT1 IFIT2 IFIT3 RSAD2 HERC5 IFI44 XAF1 BST2",
 "neuronal": "NCAM1 L1CAM GFRA1 GFRA3 NTRK2 NGFR S100B PLP1 MPZ GAP43 SOX2 STMN2 TUBB3",
 "DNA_repair": "BRCA1 BRCA2 RAD51 RAD51B XRCC2 FANCD2 FANCI ATM ATR CHEK1 CHEK2 MSH2 MLH1 PARP1 POLQ",
 "oxphos": "NDUFA4 NDUFB3 SDHA SDHB UQCRC1 UQCRC2 COX5A COX6C ATP5F1A ATP5F1B ATP5MC1 CYC1",
 "apoptosis": "BAX BAK1 CASP3 CASP7 CASP8 CASP9 BID BCL2 BCL2L1 MCL1 APAF1 PMAIP1 BBC3",
 "angiogenesis": "VEGFA VEGFB KDR FLT1 ANGPT2 TEK PECAM1 CDH5 ESM1 NOTCH4 DLL4 APLN",
 "hormone_stress": "FOS JUN JUNB EGR1 ATF3 DUSP1 KLF4 NR4A1 SOCS3 ZFP36",
}


def build(genes, ROOT):
    """Return (M, module_names, module_source) with M shape (P, K), entries in {0,1}
    signed so that +1 means 'high value of this module raises the module score'."""
    genes = list(genes)
    gi = {g: i for i, g in enumerate(genes)}
    cols, names, src = [], [], []

    def add(name, gset, source, sign=1.0):
        v = np.zeros(len(genes), np.float32)
        hit = [gi[g] for g in gset if g in gi]
        if len(hit) < 8: return
        v[hit] = sign / np.sqrt(len(hit))
        cols.append(v); names.append(name); src.append(source)

    for k, v in CURATED.items():
        add(k, v.split(), "curated")

    wg = pd.read_csv(f"{ROOT}/bulk_rna_data/outputs/images/clinical_outcome/wgcna/06_gene_info_complete.csv")
    gc = wg.columns[0]
    for mod, sub in wg.groupby("Module"):
        if mod == "grey": continue
        gl = [g for g in sub[gc].astype(str) if g in gi]
        if len(gl) < 20 or len(gl) > 2500: continue
        add(f"wgcna_{mod}", gl, "wgcna")

    for ct in ["Epithelial_cells_Tumor", "Myoepithelial_cells_Tumor", "Actively_Dividing_cells"]:
        f = f"{ROOT}/outputs/01_candidate_table/sc_pseudobulk/DE_{ct}.csv"
        de = pd.read_csv(f)
        gcol = de.columns[0]
        de = de.dropna(subset=["deseq2_padj"])
        for tag, mask in [("up_poor", (de["deseq2_padj"] < 0.05) & (de["direction"] == "Poor")),
                          ("up_good", (de["deseq2_padj"] < 0.05) & (de["direction"] == "Good"))]:
            add(f"sc_{ct}_{tag}", de.loc[mask, gcol].astype(str).tolist(), "single_cell")

    ev = pd.read_csv(f"{ROOT}/outputs/01_candidate_table/evidence/gene_evidence_wide.csv")
    for tag, d in [("poor", "Poor"), ("good", "Good")]:
        add(f"evidence_{tag}", ev.loc[ev.direction_consensus == d, "gene"].astype(str).tolist(), "evidence")
        sub = ev[(ev.direction_consensus == d) & (ev.n_sources >= 3)]
        add(f"evidence3_{tag}", sub["gene"].astype(str).tolist(), "evidence")

    M = np.stack(cols, 1)
    return M, names, src

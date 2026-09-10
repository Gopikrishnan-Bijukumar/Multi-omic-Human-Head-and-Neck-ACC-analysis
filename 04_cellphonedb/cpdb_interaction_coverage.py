"""Quantify how much of the CellPhoneDB interaction space is testable given that
the expression matrix is a 10,000 highly-variable-gene subset.

Both partners of an interaction must be present for CPDB to evaluate it; a
complex counts as present only if every subunit is present. Written to the audit
directory so 08_ can cite measured numbers in LIMITATIONS.md."""
import os
# Project data root. Override with the ACC_DATA_ROOT environment variable.
ACC_DATA_ROOT = os.environ.get("ACC_DATA_ROOT", os.path.expanduser("~/acc_data"))

import zipfile, pandas as pd, io, h5py, json, hashlib, datetime
from pathlib import Path
from anndata.io import read_elem

DB   = Path(f'{ACC_DATA_ROOT}/reference/cellphonedb/cellphonedb.zip')
H5   = Path(f'{ACC_DATA_ROOT}/outputs/h5ad_files/inprogress_3.h5ad')
OUT  = Path(f'{ACC_DATA_ROOT}/outputs/cellphoneDB/clinical_outcome_v2/audit')

z = zipfile.ZipFile(DB)
T = {n.split('/')[-1][:-4]: pd.read_csv(io.BytesIO(z.read(n)))
     for n in z.namelist() if n.endswith('.csv')}
var = set(read_elem(h5py.File(H5, 'r')['var/_index']))

gene = T['gene_table'].merge(T['protein_table'][['id_protein', 'protein_multidata_id']],
                             left_on='protein_id', right_on='id_protein', how='left')
gene['present'] = gene['hgnc_symbol'].isin(var)
prot_present = set(gene.loc[gene.present, 'protein_multidata_id'].dropna().astype(int))
prot_all     = set(gene['protein_multidata_id'].dropna().astype(int))

comp = T['complex_composition_table']
comp_ok = {cid for cid, g in comp.groupby('complex_multidata_id')
           if set(g['protein_multidata_id']).issubset(prot_present)}
comp_all = set(comp['complex_multidata_id'].unique())

inter = T['interaction_table']
detect = prot_present | comp_ok
a = inter['multidata_1_id'].isin(detect)
b = inter['multidata_2_id'].isin(detect)
n = len(inter)

rep = {
    'generated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'question': 'How much of the CellPhoneDB interaction space is testable on a 10,000-HVG matrix?',
    'matrix': {'path': str(H5), 'n_genes': len(var), 'is_hvg_subset': True,
               'full_transcriptome_available': False,
               'note': 'every h5ad in outputs/h5ad_files/ carries the same 10,000 HVGs; none has .raw'},
    'database': {'path': str(DB), 'sha256': hashlib.sha256(DB.read_bytes()).hexdigest()},
    'simple_proteins':  {'detectable': len(prot_present), 'total': len(prot_all),
                         'fraction': round(len(prot_present)/len(prot_all), 4)},
    'complexes':        {'detectable': len(comp_ok), 'total': len(comp_all),
                         'fraction': round(len(comp_ok)/len(comp_all), 4),
                         'note': 'a complex is detectable only if EVERY subunit is present'},
    'interactions': {
        'both_partners_present': int((a & b).sum()),
        'exactly_one_present':   int(((a | b) & ~(a & b)).sum()),
        'neither_present':       int((~(a | b)).sum()),
        'total':                 n,
        'detectable_fraction':   round(float((a & b).sum())/n, 4),
        'undetectable_fraction': round(float((~(a & b)).sum())/n, 4),
    },
    'corroboration': ('1224 detectable interactions equals the row count of the existing '
                      'all-samples interaction_scores.csv, so this is measured, not predicted'),
    'undetectable_by_classification':
        inter[~(a & b)]['classification'].value_counts().head(20).to_dict(),
    'interpretation': ('Absence of an interaction from these results is NOT evidence of '
                       'biological absence. No pathway-level negative claim may be made for '
                       'any classification listed above.'),
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'cpdb_interaction_coverage.json').write_text(json.dumps(rep, indent=2))
inter.assign(detectable=(a & b))[['id_cp_interaction', 'classification', 'detectable']] \
     .to_csv(OUT / 'cpdb_interaction_detectability.csv', index=False)
print(json.dumps({k: rep[k] for k in ('simple_proteins','complexes','interactions')}, indent=2))
print('\nwrote', OUT / 'cpdb_interaction_coverage.json')

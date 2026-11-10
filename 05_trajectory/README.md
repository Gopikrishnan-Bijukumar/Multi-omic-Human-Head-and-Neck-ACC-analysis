# Differentiation and trajectory analysis

Supplementary Methods, Single cell section 6. Cells labeled as malignant or carrying a tumor tag
were retained, contributing 71,984 cells from 24 single-cell patients.

| Script | What it does |
|---|---|
| `01_cytotrace2_malignant.R` | CytoTRACE2 on raw expression counts and gene counts; potency scores computed in the malignant compartment and then in all annotated cells. Pearson correlation tests the association of CytoTRACE2 scores with S and G2M phase scores |
| `02_monocle3_trajectory.R` | Monocle3 trajectory ordering on the same patient-aligned embedding, with the root selected blind to clinical outcome, plus two further independent roots to assess ordering stability |
| `03_cluster_characterization.R` | Characterizes the clusters underlying the trajectory |

### Run order

```bash
Rscript 01_cytotrace2_malignant.R      # CytoTRACE2 scoring, N_CORES=6     ~40 min
Rscript 02_monocle3_trajectory.R       # trajectory/pseudotime, cores=20   ~35-45 min
Rscript 03_cluster_characterization.R  # focal cluster + disconnected pop  ~5-10 min
```

> [!WARNING]
> **Run the stages one at a time** (peak memory is high; see the ncores note in
> `01_cytotrace2_malignant.R`). Each stage reads files the previous one writes:
>
> - after stage 1: `cytotrace_analysis/tables/cytotrace2_scores_malignant.csv`,
>   `run_config.csv`, and the cache `h5ad_files/cache/malignant_subset_counts_meta.rds`
> - after stage 2: `monocle3_analysis/tables/cds_final.rds`, `focal_cluster_id.txt`,
>   `root_annot5_composition.csv`, `annot5_class_summary.csv`
>
> Monocle3 rejects a cache holding fewer than 50,000 cells and then reads the 8 GB h5ad instead.

Stage 2 can be re-run on its own once stage 1 has written its outputs. A change to the `annot_5`
labels requires re-running stage 1, which writes them into the cache.

## `outcome_split/`

Per-outcome-cohort versions of the same analyses, plus patient-level comparisons.

| Script | What it does |
|---|---|
| `cytotrace2_cohort.R` | CytoTRACE2 per outcome cohort. Per-patient median and mean scores and the proportion of cells exceeding the cohort-wide top 10%, 25% and 33% thresholds, for samples with at least 10 malignant cells. Outcome groups compared with Wilcoxon rank-sum tests and Cliff's delta effect sizes |
| `monocle3_cohort.R` | Monocle3 per outcome cohort |
| `patient_site_distribution.R` | Patient and anatomical site distribution across the cohorts |

### Run order

Both stages take the cohort name (`good` or `poor`) as their one argument:

```bash
for c in good poor; do
  Rscript cytotrace2_cohort.R "$c"   # CytoTRACE2 on that cohort's tumor cells
  Rscript monocle3_cohort.R   "$c"   # root + trajectory + pseudotime
done
```

> [!WARNING]
> Run sequentially, not alongside the pooled pipeline above (the pooled Monocle3 run peaks near
> 48 GB). `monocle3_cohort.R` stops if the cohort cache is missing.

After each cohort, check `run_config_monocle3_<cohort>.csv` for `alignment_status`. If it reads
`FAILED_USED_PCA`, `align_cds` failed and the trajectory was built on an unaligned embedding —
every figure still renders, but that cohort's trajectory must not be interpreted.

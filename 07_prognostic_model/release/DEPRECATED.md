# Do-not-quote map

The project directory holds four generations of results. Only the files listed as **current**
in `results_dictionary.csv` are manuscript-facing. Everything below is retained for the audit
trail and must not be quoted.

> **Reading this outside the project tree.** This document was written against the working
> project directory, and most paths below (`work/...`, `outputs/...`, and the review documents
> in "Stale narrative") name **result artefacts and internal notes that are not distributed
> with this repository**. They are named so the audit trail is complete and so anyone with the
> data tree can find them — not because they ship here. The files that *do* ship are the ones
> under `release/`, plus the numbered scripts in `07_prognostic_model/`.

## Superseded model artefacts

| File | What it is | Why it is not the model |
|---|---|---|
| `outputs/04_model/final_model.json` | 30-gene elastic-net panel, 7 non-zero, apparent AUC 1.0 | Generation 1. The file self-flags `"not_evidence": true`. Its own directory records `config_path` pointing at `t2_model_build`. |
| `work/outputs/FINAL_model.json` | dense 12,700-gene Bayesian beta | Generation 2, whole-transcriptome. Superseded within hours by the compact panel. Still a live input: stages 22/27/34 read its **signs**. |
| `work/outputs/FINAL_validation.json` | DK C = 0.653 | Generation 2 validation of the dense signature, not the four-gene panel. |
| `work/outputs/validation_summary.json` | earlier generation-2 validation | Superseded by `FINAL_validation.json` the same evening. |
| `work/outputs_panel8/PANEL8_*.json` | 8-gene variant | Declared non-primary, and it **failed its own pre-specified robustness gate** (leave-one-gene-out p_max 0.333). Reported as a negative result; never promote it. |

**The model is `work/outputs_compact/COMPACT_model.json`** — DSCAM, ODC1, NCAPG, CCNB2 — restated
in `release/model_spec.json` with the frozen reference and threshold.

## Superseded numbers

| Withdrawn | Replacement | Note |
|---|---|---|
| internal discovery AUC **0.906** | 0.783 (2.5–97.5 pct 0.7475–0.8152) | The original figure was optimistic; the leakage-free procedure is the like-for-like estimate. |
| internal leakage-free **0.810** as a point estimate | **0.783 ± 0.021** over 14 draws | 0.810 is one stochastic draw and *does* lie inside the run-to-run distribution. `work/README.md`'s "0.810 is the like-for-like replacement" is one generation behind. |
| `work/outputs_metrics/INTERNAL_stability.json` and `31d_internal_stability.py` | `INTERNAL_stability_full.json` / `31e_internal_stability_full.py` | 31d concluded 0.810 was *not* reproduced; 31e, with 14 draws instead of 5, reversed that verdict. |
| `work/outputs_metrics/` (UNIFIED_metrics.json, M1–M3, S1) | `work/outputs_metrics_v2/` | v1 was generated at 09:40 on 2026-08-19, before 31e finished at 11:36, so its internal-cohort row reads 0.777 [0.761–0.795] and asserts 0.810 was not reproduced. Exactly 5 JSON keys differ; every confirmatory number is identical. |
| CCR2020 matched-null p **0.0145** (`CCR2020_results.json`) | **0.0135** (`NULLS_results.json`) | The 10,000-panel sign-structure-preserving null supersedes the earlier estimate. |
| Brayer post-hoc recovery AUC **0.955** | **0.926** | See `ACC_STUDY_CONSOLIDATED.md` revision banner. |

## Stale narrative

- `review_after_solve_delete.md` — **removed from this directory.** It was the second methods review.
  Its §"analytic transparency" quoted "mean 0.780, range 0.761–0.795, 0.810 not reproduced", which
  `31e_internal_stability_full.py` overturned 62 minutes after the review was written, so the document
  was stale on arrival on the point it pressed hardest. Its substance, and the response to every item
  in it, is preserved in `REVIEW2_RESPONSE.md`. A copy of the original survives in
  `../t4_model_build/` if the verbatim text is ever needed.
- `work/README.md` — one generation behind `ACC_STUDY_CONSOLIDATED.md` on the 0.810 question.
- `work/PRESPEC_stage1_RETROSPECTIVE.md` — retrospectively written, as its name says.

## Write hazards

- `config.yaml` sets `build_root` to **`t2_model_build`**. Running `scripts/run_all.sh` from this
  directory would write generation-1 outputs into the t2 tree. A guard has been added to
  `run_all.sh`; the paths were deliberately left alone because `outputs/01_candidate_table/` and
  `outputs/03_harmonized/` are live read-only inputs to `work/03_data.py` and must not be regenerated.
- `t4_model_build` is a byte-identical backup whose 41 `work/*.py` files still hard-code
  `ROOT = ".../t3_model_build"`. It is a snapshot, not a runnable copy: running anything from it
  reads and writes t3.

## Known duplication (left as-is deliberately)

`dwls()` exists twice, at `work/21_signature_and_gate.py:114` and `work/22_dk_composition.py:73`.
They are numerically identical. Refactoring them into a shared module mid-review would break the
"these outputs came from this code" audit trail for no numerical gain.

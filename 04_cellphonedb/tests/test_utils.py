"""Acceptance tests for clinical_outcome_v2/utils.py.

Twelve tests, one per contract that the analysis must never silently break.
All fixtures are small and in-memory: no h5ad, no CellPhoneDB, no network.
The whole file runs in seconds.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import utils


# The 19 real annot_6 labels. Hardcoding these is legitimate -- they are the
# label vocabulary the parser must survive (embedded ' - ' and '/'), not
# patient identities.
ANNOT_6_LABELS = [
    "ADC - Tumor",
    "ADMEC - Tumor",
    "Basal cells - Tumor",
    "Epithelial cells - Tumor",
    "Epithelial cells",
    "Fibroblast cells",
    "Fibroblast cells - CAF",
    "Fibroblast cells - Inflammatory",
    "Macrophages - AP/TAM",
    "Macrophages - M2",
    "CD4+ T cells",
    "CD8+ T cells",
    "B cells",
    "Memory B cells",
    "Plasma cells",
    "Mast cells",
    "Endothelial cells",
    "Mural cells",
    "Muscle cells",
]


def _write_scores_csv(path: Path, cell_pairs, rows) -> Path:
    """Write a miniature CellPhoneDB interaction_scores.csv."""
    frame = pd.DataFrame(rows)
    ordered = ["id_cp_interaction", "interacting_pair", "classification"] + list(cell_pairs)
    frame = frame.reindex(columns=ordered)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    return path


# =============================================================================
# 1. structural missingness is never turned into a zero
# =============================================================================

def test_structural_missingness_is_not_zero(tmp_path):
    """A cell type below the abundance threshold must stay unobservable."""
    sender, receiver, absent = "Fibroblast cells", "Macrophages - M2", "Mast cells"
    pairs = [
        f"{sender}|{receiver}",
        f"{sender}|{absent}",   # column exists, but 'absent' failed the threshold
    ]
    csv = _write_scores_csv(
        tmp_path / "P1" / "interaction_scores.csv",
        pairs,
        [
            {
                "id_cp_interaction": "CPI-0001",
                "interacting_pair": "A_B",
                "classification": "x",
                pairs[0]: 0.0,          # a genuinely computed ZERO
                pairs[1]: 4.2,          # score exists but the type is not retained
            }
        ],
    )
    retained = [sender, receiver]
    long_df = utils.extract_patient_endpoints(csv, retained, ANNOT_6_LABELS, patient="P1")
    out = utils.classify_missingness(long_df, {"P1": retained})

    evaluated = out.loc[out["receiver"] == receiver].iloc[0]
    structural = out.loc[out["receiver"] == absent].iloc[0]

    # a real 0.0 survives as 0.0 ...
    assert evaluated["status"] == "evaluated"
    assert evaluated["score_evaluated"] == 0.0

    # ... while a structurally missing endpoint is NaN, never 0.0
    assert structural["status"] == "structural_missing"
    assert np.isnan(structural["score_evaluated"])
    assert not out.loc[out["status"] == "structural_missing", "score_evaluated"].eq(0).any()
    assert out.loc[out["status"] == "structural_missing", "score_evaluated"].isna().all()


# =============================================================================
# 2. direction matters: A->B is a different endpoint from B->A
# =============================================================================

def test_directional_keys_distinguish_orientation(tmp_path):
    a, b = "Fibroblast cells", "Macrophages - M2"
    pairs = [f"{a}|{b}", f"{b}|{a}"]
    csv = _write_scores_csv(
        tmp_path / "P1" / "interaction_scores.csv",
        pairs,
        [{"id_cp_interaction": "CPI-0001", "interacting_pair": "A_B",
          "classification": "x", pairs[0]: 1.0, pairs[1]: 9.0}],
    )
    out = utils.classify_missingness(
        utils.extract_patient_endpoints(csv, [a, b], ANNOT_6_LABELS, patient="P1"),
        {"P1": [a, b]},
    )

    forward = out.loc[(out["sender"] == a) & (out["receiver"] == b)].iloc[0]
    reverse = out.loc[(out["sender"] == b) & (out["receiver"] == a)].iloc[0]

    assert forward["endpoint_key"] != reverse["endpoint_key"]
    assert forward["endpoint_key"] == utils.make_endpoint_key("CPI-0001", a, b)
    assert forward["score_evaluated"] == 1.0
    assert reverse["score_evaluated"] == 9.0
    assert out["endpoint_key"].is_unique

    # and the parser keeps orientation too
    assert utils.parse_cellpair_columns(pairs, ANNOT_6_LABELS) == [(a, b), (b, a)]


# =============================================================================
# 3. duplicate endpoint rows fail loudly (flagged invalid, never averaged away)
# =============================================================================

def test_duplicate_endpoint_rows_fail_loudly(tmp_path):
    a, b = "Fibroblast cells", "Macrophages - M2"
    pairs = [f"{a}|{b}"]
    csv = _write_scores_csv(
        tmp_path / "P1" / "interaction_scores.csv",
        pairs,
        [
            {"id_cp_interaction": "CPI-0001", "interacting_pair": "A_B",
             "classification": "x", pairs[0]: 1.0},
            {"id_cp_interaction": "CPI-0001", "interacting_pair": "A_B",
             "classification": "x", pairs[0]: 7.0},   # duplicate id
        ],
    )
    out = utils.classify_missingness(
        utils.extract_patient_endpoints(csv, [a, b], ANNOT_6_LABELS, patient="P1"),
        {"P1": [a, b]},
    )
    dupes = out.loc[out["endpoint_key"] == utils.make_endpoint_key("CPI-0001", a, b)]
    assert len(dupes) == 2
    assert (dupes["status"] == "invalid").all()
    assert dupes["score_evaluated"].isna().all()
    assert not dupes["tested"].any()

    # a duplicate is never quietly collapsed into one evaluated endpoint
    assert (out["status"] == "evaluated").sum() == 0


# =============================================================================
# 4. Mann-Whitney effect sign, hand-checkable
# =============================================================================

def test_mannwhitney_effect_sign_is_correct():
    good = [5.0, 6.0, 7.0]
    poor = [1.0, 2.0, 3.0]

    U, p, method = utils.mannwhitney_hybrid(good, poor)
    assert method == "exact"
    assert U == 9.0                                   # every good > every poor: 3 x 3
    assert utils.rank_biserial(U, 3, 3) == pytest.approx(1.0)

    U_rev, p_rev, _ = utils.mannwhitney_hybrid(poor, good)
    assert U_rev == 0.0
    assert utils.rank_biserial(U_rev, 3, 3) == pytest.approx(-1.0)
    assert p == pytest.approx(p_rev)                  # two-sided p is symmetric

    # partial overlap, still hand-checkable: good=[2,4,6] vs poor=[1,3,5]
    U_mix, _, _ = utils.mannwhitney_hybrid([2.0, 4.0, 6.0], [1.0, 3.0, 5.0])
    assert U_mix == 6.0                               # 1 + 2 + 3 wins
    assert utils.rank_biserial(U_mix, 3, 3) == pytest.approx(2 * 6 / 9 - 1)

    # sign convention: positive means larger in Good
    assert utils.rank_biserial(*[utils.mannwhitney_hybrid(good, poor)[0], 3, 3]) > 0


# =============================================================================
# 5. BH is applied to all and only the tested endpoints
# =============================================================================

def test_bh_applied_to_tested_endpoints_only():
    frame = pd.DataFrame(
        {
            "p": [0.001, 0.02, 0.30, np.nan, 0.75, np.nan],
            "tested": [True, True, True, False, True, False],
        }
    )
    # contract: untested endpoints enter bh() as NaN
    frame["q"] = utils.bh(frame["p"].to_numpy())

    assert frame.loc[~frame["tested"], "q"].isna().all()
    assert frame.loc[frame["tested"], "q"].notna().all()

    # family size is the number of TESTED endpoints (4), not the frame length (6)
    tested_p = np.sort(frame.loc[frame["tested"], "p"].to_numpy())
    m = tested_p.size
    assert m == 4
    expected = np.minimum.accumulate((tested_p * m / np.arange(1, m + 1))[::-1])[::-1]
    got = np.sort(frame.loc[frame["tested"], "q"].to_numpy())
    assert got == pytest.approx(expected)

    # sanity: had the NaNs been counted, m would be 6 and q would be larger
    assert frame.loc[0, "q"] == pytest.approx(0.001 * 4 / 1)
    assert frame.loc[0, "q"] < 0.001 * 6

    # untested endpoints must not shift the tested ones
    only_tested = utils.bh(frame.loc[frame["tested"], "p"].to_numpy())
    assert np.sort(only_tested) == pytest.approx(expected)


# =============================================================================
# 6. the bootstrap resamples patients within outcome, never cells
# =============================================================================

def test_bootstrap_resamples_patients_within_outcome():
    n_good, n_poor = 13, 11
    rng = np.random.default_rng(0)
    # every Good patient == 1.0, every Poor patient == 0.0. Any resample that
    # stays within its own outcome group can only ever reproduce md = 1.0 and
    # rb = 1.0; mixing patients across groups, or resampling at the cell level,
    # would widen the interval.
    values = np.concatenate(
        [np.ones((5, n_good)), np.zeros((5, n_poor))], axis=1
    )
    good_idx = np.arange(n_good)
    poor_idx = np.arange(n_good, n_good + n_poor)

    out = utils.bootstrap_ci(values, good_idx, poor_idx, n_resamples=500, seed=7)

    assert out["good_resample_shape"] == (500, n_good)   # patients, not cells
    assert out["poor_resample_shape"] == (500, n_poor)
    assert np.allclose(out["median_difference"], 1.0)
    assert np.allclose(out["median_difference_ci_low"], 1.0)
    assert np.allclose(out["median_difference_ci_high"], 1.0)
    assert np.allclose(out["rank_biserial"], 1.0)
    assert np.allclose(out["rank_biserial_ci_low"], 1.0)

    # a heterogeneous endpoint DOES get a non-degenerate interval, so the
    # zero-width result above is a property of the design, not of a no-op
    het = rng.normal(size=(3, n_good + n_poor))
    het_out = utils.bootstrap_ci(het, good_idx, poor_idx, n_resamples=500, seed=7)
    assert np.all(het_out["median_difference_ci_high"] > het_out["median_difference_ci_low"])

    # the same seed reproduces the same world; a different seed does not have to
    again = utils.bootstrap_ci(het, good_idx, poor_idx, n_resamples=500, seed=7)
    assert np.allclose(het_out["median_difference_ci_low"], again["median_difference_ci_low"])

    # NaN (structurally missing patient) is tolerated, not silently zeroed
    with_nan = values.copy()
    with_nan[0, 0] = np.nan
    nan_out = utils.bootstrap_ci(with_nan, good_idx, poor_idx, n_resamples=200, seed=3)
    assert nan_out["n_good_evaluable"][0] == n_good - 1
    assert np.isclose(nan_out["median_difference"][0], 1.0)


# =============================================================================
# 7. run_is_complete rejects missing outputs and stale configs
# =============================================================================

def test_run_is_complete_rejects_missing_file_or_stale_hash(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    required = ["a.csv", "b.csv"]
    for name in required:
        (run_dir / name).write_text("x")

    cfg_hash = utils.config_hash({"a": 1, "b": [2, 3]})
    assert not utils.run_is_complete(run_dir, required, cfg_hash)   # no _SUCCESS yet

    utils.write_success_json(run_dir, config_hash=cfg_hash, note="ok")
    assert utils.run_is_complete(run_dir, required, cfg_hash)

    # one missing output
    (run_dir / "b.csv").unlink()
    assert not utils.run_is_complete(run_dir, required, cfg_hash)
    (run_dir / "b.csv").write_text("x")
    assert utils.run_is_complete(run_dir, required, cfg_hash)

    # mismatched config hash
    other_hash = utils.config_hash({"a": 1, "b": [2, 4]})
    assert other_hash != cfg_hash
    assert not utils.run_is_complete(run_dir, required, other_hash)

    # config_hash ignores key order but tracks values
    assert utils.config_hash({"b": [2, 3], "a": 1}) == cfg_hash

    # and _SUCCESS.json is real JSON carrying the hash
    record = json.loads((run_dir / utils.SUCCESS_FILENAME).read_text())
    assert record["config_hash"] == cfg_hash
    assert "completed_utc" in record


# =============================================================================
# 8. one patient may not carry two outcomes
# =============================================================================

def test_derive_patient_outcomes_rejects_conflicting_outcome():
    good_obs = pd.DataFrame(
        {
            "batch": ["P1", "P1", "P2", "P2", "P3"],
            "clinical_outcome": ["Good", "Good", "Poor", "Poor", "Good"],
        },
        index=[f"c{i}" for i in range(5)],
    )
    outcomes = utils.derive_patient_outcomes(good_obs, "batch", "clinical_outcome")
    assert outcomes.to_dict() == {"P1": "Good", "P2": "Poor", "P3": "Good"}
    assert list(outcomes.index) == ["P1", "P2", "P3"]      # sorted, no hardcoded ids

    bad_obs = good_obs.copy()
    bad_obs.loc["c1", "clinical_outcome"] = "Poor"          # P1 now Good AND Poor
    with pytest.raises(ValueError, match="more than one"):
        utils.derive_patient_outcomes(bad_obs, "batch", "clinical_outcome")


# =============================================================================
# 9. count validation never densifies
# =============================================================================

def test_count_validation_on_sparse_without_densifying():
    class NoDensify(sp.csr_matrix):
        """A CSR matrix that explodes if anything tries to densify it."""

        def toarray(self, *args, **kwargs):          # pragma: no cover - must not run
            raise AssertionError("densified via toarray()")

        def todense(self, *args, **kwargs):          # pragma: no cover - must not run
            raise AssertionError("densified via todense()")

        def __array__(self, *args, **kwargs):        # pragma: no cover - must not run
            raise AssertionError("densified via __array__()")

    dense = np.zeros((200, 300), dtype=np.float32)
    dense[0, 0] = 3
    dense[5, 7] = 11
    counts = NoDensify(sp.csr_matrix(dense))

    summary = utils.validate_counts_matrix(counts, name="counts")
    assert summary["shape"] == (200, 300)
    assert summary["n_stored"] == 2            # only the stored non-zeros were touched
    assert summary["sum"] == 14.0
    assert summary["densified"] is False

    # non-integral (z-scored) data is rejected, still without densifying
    zscored_dense = np.zeros((10, 10), dtype=np.float32)
    zscored_dense[1, 1] = 0.5
    zscored = NoDensify(sp.csr_matrix(zscored_dense))
    with pytest.raises(ValueError, match="not integral"):
        utils.validate_counts_matrix(zscored, name="counts")

    negative_dense = np.zeros((10, 10), dtype=np.float32)
    negative_dense[2, 2] = -3
    with pytest.raises(ValueError, match="negative"):
        utils.validate_counts_matrix(NoDensify(sp.csr_matrix(negative_dense)))


# =============================================================================
# 10. the parser survives all 19 real annot_6 labels
# =============================================================================

def test_parse_cellpair_columns_with_all_real_labels():
    assert len(ANNOT_6_LABELS) == 19
    assert not any("|" in label for label in ANNOT_6_LABELS)

    cell_pair_cols = [f"{s}|{r}" for s in ANNOT_6_LABELS for r in ANNOT_6_LABELS]
    columns = ["id_cp_interaction", "interacting_pair", "classification",
               "is_integrin", "directionality"] + cell_pair_cols

    pairs = utils.parse_cellpair_columns(columns, ANNOT_6_LABELS)

    assert len(pairs) == 19 * 19 == 361                    # metadata columns skipped
    assert pairs == [(s, r) for s in ANNOT_6_LABELS for r in ANNOT_6_LABELS]
    # labels containing ' - ' and '/' round-trip intact
    assert ("Macrophages - AP/TAM", "Fibroblast cells - Inflammatory") in pairs
    assert ("Epithelial cells", "Epithelial cells - Tumor") in pairs
    assert ("Epithelial cells - Tumor", "Epithelial cells") in pairs

    # an unknown half raises loudly rather than being dropped
    with pytest.raises(ValueError, match="unknown cell type"):
        utils.parse_cellpair_columns(["Fibroblast cells|Not A Cell Type"], ANNOT_6_LABELS)
    # so does a malformed column with two separators
    with pytest.raises(ValueError, match="expected 2"):
        utils.parse_cellpair_columns(["B cells|B cells|B cells"], ANNOT_6_LABELS)


# =============================================================================
# 11. safe_unlink refuses to escape its root
# =============================================================================

def test_safe_unlink_refuses_paths_outside_root(tmp_path):
    root = tmp_path / "outputs"
    root.mkdir()
    inside = root / "figure.png"
    inside.write_text("keep-until-asked")
    outside = tmp_path / "precious.h5ad"
    outside.write_text("do not delete")

    # escape by absolute path
    with pytest.raises(ValueError, match="outside the allowed root"):
        utils.safe_unlink(outside, root)
    assert outside.exists()

    # escape by traversal
    with pytest.raises(ValueError, match="outside the allowed root"):
        utils.safe_unlink(root / ".." / "precious.h5ad", root)
    assert outside.exists()

    # escape by symlink
    link = root / "link.h5ad"
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="outside the allowed root"):
        utils.safe_unlink(link, root)
    assert outside.exists()

    # the root itself is protected, and directories are refused
    with pytest.raises(ValueError, match="allowed root itself"):
        utils.safe_unlink(root, root)
    subdir = root / "sub"
    subdir.mkdir()
    with pytest.raises(IsADirectoryError):
        utils.safe_unlink(subdir, root)

    # a legitimate file inside the root IS removed, and a missing file is a no-op
    utils.safe_unlink(inside, root)
    assert not inside.exists()
    utils.safe_unlink(inside, root)


# =============================================================================
# 12. the exact/asymptotic hybrid, and why it matters at 13 vs 11
# =============================================================================

def test_mannwhitney_hybrid_exact_when_tiefree_asymptotic_when_tied():
    n_good, n_poor = 13, 11

    # perfect separation, no ties anywhere -> exact null
    good = np.arange(n_good, dtype=float) + 100.0
    poor = np.arange(n_poor, dtype=float)
    U, p, method = utils.mannwhitney_hybrid(good, poor)
    assert method == "exact"
    assert U == float(n_good * n_poor)
    assert utils.rank_biserial(U, n_good, n_poor) == pytest.approx(1.0)
    assert p == pytest.approx(8.012358e-07, rel=1e-5)

    # the asymptotic null cannot resolve below ~3.90e-05 at this n ...
    from scipy.stats import mannwhitneyu

    p_asymptotic = mannwhitneyu(good, poor, alternative="two-sided", method="asymptotic").pvalue
    assert p_asymptotic == pytest.approx(3.897139e-05, rel=1e-5)
    # ... and under BH over ~20,000 endpoints that floor can never reach q<0.05,
    # while the exact value can. This is the whole reason for the hybrid.
    family = 20_000
    assert p_asymptotic * family > 0.05
    assert p * family < 0.05

    # a single tie flips the method to asymptotic (exact is invalid under ties)
    tied_good = good.copy()
    tied_good[0] = poor[0]
    U_t, p_t, method_t = utils.mannwhitney_hybrid(tied_good, poor)
    assert method_t == "asymptotic"
    assert np.isfinite(p_t)

    # a tie WITHIN one group also counts
    within_tie = good.copy()
    within_tie[1] = within_tie[0]
    assert utils.mannwhitney_hybrid(within_tie, poor)[2] == "asymptotic"

    # empty group after dropping non-finite values -> declared insufficient
    U_n, p_n, method_n = utils.mannwhitney_hybrid([np.nan, np.nan], poor)
    assert method_n == "insufficient"
    assert np.isnan(U_n) and np.isnan(p_n)

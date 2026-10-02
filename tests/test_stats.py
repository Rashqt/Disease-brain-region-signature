"""All data in this file is SYNTHETIC test data."""
import numpy as np
import pandas as pd

from dbrs.stats.null import bh_adjust, draw_matched_sets, empirical_p, run_signature


def synthetic_matrix(n_regions=4, n_genes=200, seed=0):
    rng = np.random.default_rng(seed)
    z = pd.DataFrame(rng.normal(size=(n_regions, n_genes)),
                     index=range(n_regions), columns=[f"G{i}" for i in range(n_genes)])
    level = pd.Series(rng.normal(size=n_genes), index=z.columns)
    return z, level


def test_bh_matches_hand_calculation():
    adj = bh_adjust([0.01, 0.04, 0.03, 0.2])
    assert np.allclose(adj, [0.04, 0.0533333, 0.0533333, 0.2])


def test_bh_matches_statsmodels_if_installed():
    sm = __import__("pytest").importorskip("statsmodels.stats.multitest")
    p = np.random.default_rng(1).uniform(size=30)
    assert np.allclose(bh_adjust(p), sm.multipletests(p, method="fdr_bh")[1])


def test_empirical_p_never_zero_and_correct():
    null = np.arange(100, dtype=float)[None, :]
    assert empirical_p(np.array([1000.0]), null)[0] == 1 / 101
    assert empirical_p(np.array([50.0]), null)[0] == (1 + 50) / 101


def test_matched_sets_exclude_disease_genes_and_respect_bins():
    z, level = synthetic_matrix()
    idx = np.arange(10)
    sets = draw_matched_sets(level, idx, 50, 5, np.random.default_rng(3))
    assert not np.isin(sets, idx).any()
    assert all(len(set(row)) == len(row) for row in sets)  # no duplicates within a set
    from dbrs.stats.null import level_bins
    bins = level_bins(level, 5)
    assert (bins[sets] == bins[idx][None, :]).all()


def test_signature_is_reproducible_for_a_seed():
    z, level = synthetic_matrix()
    genes = [f"G{i}" for i in range(12)]
    a, _ = run_signature(z, level, genes, 200, 100, 5, seed=42)
    b, _ = run_signature(z, level, genes, 200, 100, 5, seed=42)
    c, _ = run_signature(z, level, genes, 200, 100, 5, seed=43)
    assert a.equals(b)
    assert not a["null_mean"].equals(c["null_mean"])


def test_planted_signal_is_detected_and_noise_is_not():
    z, level = synthetic_matrix(n_genes=400)
    genes = [f"G{i}" for i in range(20)]
    z.loc[0, genes] += 2.0                                  # planted synthetic enrichment in region 0
    res, _ = run_signature(z, level, genes, 1000, 100, 5, seed=1)
    assert res.loc[res.region_id == 0, "p_empirical"].iloc[0] < 0.01
    assert (res.loc[res.region_id != 0, "p_adj_bh"] > 0.05).all()

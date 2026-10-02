"""Matched-null model for regional gene-set scores.

Observed score for a region = mean regional z-score over the disease genes.
Null = the same score for random gene sets of equal size, drawn without replacement within a set
from genes NOT in the disease set, matched gene-by-gene on expression-level bin.
Matching on gene length / other confounders is NOT implemented yet (see README limitations).
"""
import numpy as np
import pandas as pd


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    scaled = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def level_bins(level, n_bins):
    ranks = level.rank(method="first")
    return pd.qcut(ranks, n_bins, labels=False).to_numpy()


def draw_matched_sets(level, gene_idx, n_perm, n_bins, rng):
    """Positions (into level.index) of matched random genes: array (n_perm, len(gene_idx))."""
    bins = level_bins(level, n_bins)
    in_set = np.zeros(len(level), dtype=bool)
    in_set[gene_idx] = True
    gene_idx = np.asarray(gene_idx)
    out = np.empty((n_perm, len(gene_idx)), dtype=np.int64)
    for b in np.unique(bins[gene_idx]):
        slots = np.where(bins[gene_idx] == b)[0]
        pool = np.where((bins == b) & ~in_set)[0]
        k = len(slots)
        if len(pool) < k:
            raise ValueError(f"Expression bin {b} has {len(pool)} candidate genes but {k} are needed.")
        picks = rng.random((n_perm, len(pool))).argpartition(k - 1, axis=1)[:, :k]
        out[:, slots] = pool[picks]
    return out


def score_sets(z, sets, chunk=500):
    """z: regions x genes ndarray. sets: (n_sets, n_genes) positions. Returns (regions, n_sets)."""
    out = np.empty((z.shape[0], sets.shape[0]))
    for start in range(0, sets.shape[0], chunk):
        block = sets[start:start + chunk]
        out[:, start:start + chunk] = z[:, block].mean(axis=2)
    return out


def empirical_p(observed, null, alternative="greater"):
    """observed: (R,), null: (R, N). Add-one empirical p so it is never exactly 0."""
    n = null.shape[1]
    if alternative == "greater":
        hits = (null >= observed[:, None]).sum(axis=1)
    elif alternative == "two-sided":
        centered = np.abs(null - null.mean(axis=1, keepdims=True))
        hits = (centered >= np.abs(observed - null.mean(axis=1))[:, None]).sum(axis=1)
    else:
        raise ValueError("alternative must be 'greater' or 'two-sided'")
    return (1 + hits) / (n + 1)


def run_signature(z_df, level, gene_set, n_perm, n_bootstrap, n_bins, seed, alternative="greater"):
    """Return (results DataFrame, null score matrix). Genes in gene_set must exist in z_df columns."""
    genes = [g for g in gene_set if g in z_df.columns]
    if len(genes) < 2:
        raise ValueError("Need at least 2 disease genes present in the atlas matrix.")
    level = level.reindex(z_df.columns)
    z = z_df.to_numpy()
    gene_idx = np.array([z_df.columns.get_loc(g) for g in genes])

    seeds = np.random.SeedSequence(seed).spawn(2)
    rng_null, rng_boot = (np.random.default_rng(s) for s in seeds)

    observed = z[:, gene_idx].mean(axis=1)
    null = score_sets(z, draw_matched_sets(level, gene_idx, n_perm, n_bins, rng_null))
    boot_sets = gene_idx[rng_boot.integers(0, len(gene_idx), size=(n_bootstrap, len(gene_idx)))]
    boot = score_sets(z, boot_sets)

    null_mean, null_sd = null.mean(axis=1), null.std(axis=1, ddof=1)
    p = empirical_p(observed, null, alternative)
    results = pd.DataFrame({
        "region_id": z_df.index,
        "n_genes": len(genes),
        "observed": observed,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "effect_diff": observed - null_mean,
        "effect_z": np.where(null_sd > 0, (observed - null_mean) / np.where(null_sd > 0, null_sd, 1), np.nan),
        "boot_ci_low": np.percentile(boot, 2.5, axis=1),
        "boot_ci_high": np.percentile(boot, 97.5, axis=1),
        "p_empirical": p,
        "p_mc_se": np.sqrt(p * (1 - p) / n_perm),
        "p_adj_bh": bh_adjust(p),
    })
    return results, null

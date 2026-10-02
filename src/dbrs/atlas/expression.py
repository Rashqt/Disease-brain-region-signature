"""Probe filtering, probe->gene collapsing and regional expression matrices."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def present_probes(pacalls, min_fraction):
    """Probes with a present call in >= min_fraction of all pooled samples."""
    pooled = pd.concat(pacalls, axis=1)
    return pooled.index[pooled.mean(axis=1) >= min_fraction]


def region_means(expr, region_ids, probe_ids, min_samples=1):
    """Mean expression per region for the given probes. Regions with too few samples are dropped."""
    sub = expr.loc[probe_ids]
    regions = pd.Series(np.asarray(region_ids))
    counts = regions.value_counts()
    keep_regions = counts.index[counts >= min_samples]
    grouped = sub.T.groupby(regions.values).mean().T
    return grouped[[r for r in grouped.columns if r in set(keep_regions)]]


def differential_stability(means_by_donor, min_shared=5):
    """Mean pairwise Spearman correlation of a probe's regional profile across donors."""
    probes = means_by_donor[0].index
    total = pd.Series(0.0, index=probes)
    n_pairs = 0
    for i in range(len(means_by_donor)):
        for j in range(i + 1, len(means_by_donor)):
            a, b = means_by_donor[i], means_by_donor[j]
            shared = [c for c in a.columns if c in set(b.columns)]
            if len(shared) < min_shared:
                continue
            ra = rankdata(a[shared].to_numpy(), axis=1)
            rb = rankdata(b[shared].to_numpy(), axis=1)
            ra = ra - ra.mean(axis=1, keepdims=True)
            rb = rb - rb.mean(axis=1, keepdims=True)
            denom = np.sqrt((ra ** 2).sum(axis=1) * (rb ** 2).sum(axis=1))
            with np.errstate(invalid="ignore", divide="ignore"):
                corr = (ra * rb).sum(axis=1) / denom
            total += np.nan_to_num(corr, nan=-1.0)
            n_pairs += 1
    if n_pairs == 0:
        raise ValueError("No donor pair shares enough regions to compute differential stability.")
    return total / n_pairs


def select_probe_per_gene(probe_genes, stability):
    """Keep the probe with the highest differential stability for each gene symbol."""
    df = pd.DataFrame({"gene": probe_genes, "ds": stability.reindex(probe_genes.index)})
    df = df.sort_values(["gene", "ds"], ascending=[True, False], kind="mergesort")
    return df.drop_duplicates("gene", keep="first")


def zscore_rows(df):
    sd = df.std(axis=1, ddof=1)
    return df.sub(df.mean(axis=1), axis=0).div(sd.where(sd > 0), axis=0)


def regional_gene_matrix(gene_means_by_donor, min_donors):
    """Combine per-donor gene x region means into region x gene z-scores plus a per-gene level.

    z-scoring is per donor, per gene, across that donor's regions (each region weighs equally,
    so densely sampled cortex does not dominate). Region z-scores are then averaged over donors.
    """
    zs = [zscore_rows(m) for m in gene_means_by_donor]
    donor_counts = pd.concat([pd.Series(1, index=m.columns) for m in gene_means_by_donor], axis=1).sum(axis=1)
    kept = donor_counts.index[donor_counts >= min_donors]
    z = pd.concat(zs, axis=0).groupby(level=0).mean()
    z = z[[c for c in z.columns if c in set(kept)]]
    level = pd.concat(gene_means_by_donor, axis=0).groupby(level=0).mean().mean(axis=1)
    z = z.dropna(axis=0, how="any")
    z_regions = z.T.sort_index()
    return z_regions, level.reindex(z_regions.columns), donor_counts.reindex(z_regions.index)

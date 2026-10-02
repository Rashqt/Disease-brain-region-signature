import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import data_path


def _label(row):
    return f"{row['acronym']} ({int(row['n_samples'])})" if pd.notna(row.get("acronym")) else str(row["region_id"])


def plot_signature(results, alpha, path):
    df = results.sort_values("effect_z")
    fig, ax = plt.subplots(figsize=(7, max(3, 0.28 * len(df))))
    colors = np.where(df["p_adj_bh"] < alpha, "#c0392b", "#7f8c8d")
    ax.scatter(df["effect_z"], range(len(df)), c=colors, s=18)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([_label(r) for _, r in df.iterrows()], fontsize=7)
    ax.set_xlabel("(observed - null mean) / null SD")
    ax.set_title(f"Regional signature; red = BH-adjusted p < {alpha}\nlabel: acronym (n samples)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_null(results, null, region_ids, path, top=6):
    index = {rid: i for i, rid in enumerate(region_ids)}
    top_rows = results.sort_values("effect_z", ascending=False).head(top)
    fig, axes = plt.subplots(2, 3, figsize=(10, 5), squeeze=False)
    for ax, (_, row) in zip(axes.ravel(), top_rows.iterrows()):
        ax.hist(null[index[row["region_id"]]], bins=50, color="#bdc3c7")
        ax.axvline(row["observed"], color="#c0392b")
        ax.set_title(f"{_label(row)}  p={row['p_empirical']:.4f}", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Observed (red) vs matched-null; top {top} regions by effect size (not by prior expectation)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_heatmap(z, genes, results, path):
    order = results.sort_values("effect_z", ascending=False)["region_id"].tolist()
    sub = z.loc[order, genes]
    sub = sub[sub.mean(axis=0).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(max(6, 0.18 * sub.shape[1]), max(4, 0.2 * sub.shape[0])))
    im = ax.imshow(sub.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(sub.shape[1]))
    ax.set_xticklabels(sub.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(sub.shape[0]))
    ax.set_yticklabels(sub.index, fontsize=6)
    fig.colorbar(im, ax=ax, label="regional z-score")
    ax.set_title("Disease-gene regional expression (rows ordered by effect size; row labels = region_id)", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_brain_scatter(samples, results, path):
    """Real AHBA sample coordinates (Allen-provided MNI), coloured by region effect. Regions without a result are omitted."""
    merged = samples.merge(results[["region_id", "effect_z"]], on="region_id", how="inner")
    lim = np.nanpercentile(np.abs(merged["effect_z"]), 98) or 1.0
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, (a, b, title) in zip(axes, [("mni_y", "mni_z", "sagittal (y, z)"), ("mni_x", "mni_y", "axial (x, y)"),
                                         ("mni_x", "mni_z", "coronal (x, z)")]):
        sc = ax.scatter(merged[a], merged[b], c=merged["effect_z"], cmap="RdBu_r", vmin=-lim, vmax=lim, s=5)
        ax.set_title(title, fontsize=9)
        ax.set_aspect("equal")
    fig.colorbar(sc, ax=axes, label="region effect (null z)", shrink=0.8)
    fig.suptitle("AHBA tissue samples coloured by their region's effect. Allen-provided MNI coordinates; not a surface/volume map.", fontsize=8)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def make_all(cfg, results, null, z, usable):
    out = data_path(cfg, "results", "figures", ".keep").parent
    alpha = cfg["stats"]["alpha"]
    region_ids = np.load(data_path(cfg, "results", "null_scores.npz"))["region_id"]
    plot_signature(results, alpha, out / "signature_by_region.png")
    plot_null(results, null, region_ids, out / "observed_vs_null.png")
    plot_heatmap(z, usable, results, out / "gene_region_heatmap.png")
    samples = pd.read_csv(data_path(cfg, "data", "processed", "sample_regions.csv.gz"))
    plot_brain_scatter(samples, results, out / "brain_samples_scatter.png")

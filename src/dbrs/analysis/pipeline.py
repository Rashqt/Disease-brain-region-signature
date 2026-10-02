import numpy as np
import pandas as pd

from ..atlas import expression as ex
from ..atlas import regions as rg
from ..atlas.load import load_all_donors
from ..config import data_path
from ..gwas import client, genes, parse
from ..provenance import record


def gwas_stage(cfg, refresh=False):
    recs, meta = client.fetch_associations(cfg, refresh=refresh)
    assoc, fields = parse.parse_associations(recs)
    assoc.to_csv(data_path(cfg, "data", "processed", "gwas_associations.csv"), index=False)
    table = genes.build_gene_table(assoc, cfg["gwas"]["p_threshold"], tuple(cfg["gwas"]["gene_sources"]))
    table.to_csv(data_path(cfg, "data", "processed", "disease_genes.csv"), index=False)
    record(cfg, "gwas", **meta, field_mapping=fields, p_threshold=cfg["gwas"]["p_threshold"],
           gene_sources=cfg["gwas"]["gene_sources"], n_rows=len(assoc), n_genes=len(table))
    return assoc, table


def atlas_stage(cfg):
    a = cfg["atlas"]
    donors = load_all_donors(cfg)
    ontology = donors[0].ontology
    pa = ex.present_probes([d.pacall for d in donors], a["min_present_fraction"])
    probe_info = donors[0].probes
    keep = probe_info.index[probe_info.index.isin(pa) & probe_info["entrez_id"].notna()
                            & probe_info["gene_symbol"].notna()]

    region_tables, probe_means, sample_regions = [], [], []
    for d in donors:
        ids, table = rg.build_region_map(d.samples, ontology, a["region_depth"])
        sample_regions.append(pd.DataFrame({"donor": d.donor_id, "region_id": ids,
                                            "mni_x": d.samples["mni_x"], "mni_y": d.samples["mni_y"],
                                            "mni_z": d.samples["mni_z"],
                                            "structure_acronym": d.samples["structure_acronym"]}))
        region_tables.append(table)
        probe_means.append(ex.region_means(d.expr, ids, keep, a["min_samples_per_donor_region"]))

    stability = ex.differential_stability(probe_means, a["min_shared_regions_for_ds"])
    symbols = probe_info.loc[keep, "gene_symbol"].astype(str).str.upper()
    chosen = ex.select_probe_per_gene(symbols, stability)
    gene_means = []
    for m in probe_means:
        gm = m.loc[chosen.index]
        gm.index = chosen["gene"].values
        gene_means.append(gm)

    z, level, donor_counts = ex.regional_gene_matrix(gene_means, a["min_donors_per_region"])
    region_info = pd.concat(region_tables).drop_duplicates("region_id").set_index("region_id")
    sample_df = pd.concat(sample_regions, ignore_index=True)
    n_samples = sample_df.groupby("region_id").size()
    meta = region_info.reindex(z.index).assign(n_donors=donor_counts, n_samples=n_samples.reindex(z.index))
    meta.to_csv(data_path(cfg, "data", "processed", "regions.csv"))
    z.to_csv(data_path(cfg, "data", "processed", "regional_expression_z.csv.gz"))
    level.rename("level").to_csv(data_path(cfg, "data", "processed", "gene_level.csv"))
    sample_df.to_csv(data_path(cfg, "data", "processed", "sample_regions.csv.gz"), index=False)
    record(cfg, "atlas", donors=a["donors"], region_depth=a["region_depth"], n_probes_kept=len(keep),
           n_genes=z.shape[1], n_regions=z.shape[0], min_present_fraction=a["min_present_fraction"],
           min_donors_per_region=a["min_donors_per_region"],
           note="AHBA copyright Allen Institute; raw files are not redistributed here.")
    return z, level, meta


def map_gene_set(disease_genes, atlas_genes):
    """Exact HGNC-symbol match. Unmatched genes are reported, never silently remapped or dropped."""
    atlas = set(atlas_genes)
    report = disease_genes.copy()
    report["in_atlas"] = report["gene"].isin(atlas)
    return report


def signature_stage(cfg, z=None, level=None):
    from ..stats.null import run_signature

    s = cfg["stats"]
    if z is None:
        z = pd.read_csv(data_path(cfg, "data", "processed", "regional_expression_z.csv.gz"), index_col=0)
        z.index = z.index.astype(int)
        level = pd.read_csv(data_path(cfg, "data", "processed", "gene_level.csv"), index_col=0)["level"]
    disease = pd.read_csv(data_path(cfg, "data", "processed", "disease_genes.csv"))
    report = map_gene_set(disease, z.columns)
    report.to_csv(data_path(cfg, "data", "processed", "gene_set_mapping_report.csv"), index=False)
    usable = report.loc[report["in_atlas"], "gene"].tolist()

    results, null = run_signature(z, level, usable, s["n_permutations"], s["n_bootstrap"],
                                  s["n_level_bins"], s["seed"], s["alternative"])
    meta = pd.read_csv(data_path(cfg, "data", "processed", "regions.csv"), index_col=0)
    results = results.merge(meta.reset_index(), on="region_id", how="left")
    results = results.sort_values("p_empirical").reset_index(drop=True)
    results.to_csv(data_path(cfg, "results", "signature_results.csv"), index=False)
    np.savez_compressed(data_path(cfg, "results", "null_scores.npz"), null=null, region_id=z.index.to_numpy())
    record(cfg, "signature", n_disease_genes=len(disease), n_genes_in_atlas=len(usable),
           n_genes_not_in_atlas=int((~report["in_atlas"]).sum()), **{f"stats_{k}": v for k, v in s.items()})
    return results, null, z, usable

import tomllib

import pandas as pd

from ..config import data_path


def load_known_regions(cfg):
    path = data_path(cfg, "data", "metadata", "known_vulnerable_regions.toml")
    with open(path, "rb") as fh:
        known = tomllib.load(fh)
    if not known.get("source", "").strip() or not known.get("region"):
        raise SystemExit("known_vulnerable_regions.toml needs a cited `source` and at least one [[region]]. "
                         "Fill it from a real source before comparing (and not after peeking at results).")
    return known


def compare(cfg, results):
    """Look up each independently documented region in the computed signature. No tuning, no new test."""
    known = load_known_regions(cfg)
    samples = pd.read_csv(data_path(cfg, "data", "processed", "sample_regions.csv.gz"))
    rows = []
    for item in known["region"]:
        acronyms = set(item.get("atlas_acronyms", []))
        region_ids = sorted(samples.loc[samples["structure_acronym"].isin(acronyms), "region_id"].unique())
        hit = results[results["region_id"].isin(region_ids)]
        if hit.empty:
            rows.append({"known_region": item["name"], "atlas_region": None, "status": "not covered by atlas regions used",
                         "note": item.get("note", "")})
            continue
        for _, r in hit.iterrows():
            rows.append({"known_region": item["name"], "atlas_region": r["name"], "effect_z": r["effect_z"],
                         "p_empirical": r["p_empirical"], "p_adj_bh": r["p_adj_bh"],
                         "status": "covered", "note": item.get("note", "")})
    out = pd.DataFrame(rows)
    out.to_csv(data_path(cfg, "results", "known_region_comparison.csv"), index=False)
    return out, known["source"]

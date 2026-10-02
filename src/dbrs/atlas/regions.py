"""Sample -> region mapping from the Allen structure ontology.

A sample's structure_id_path is a chain of ancestors. A coarse region is the ancestor at
`depth` (0-based position in the path); samples from structures shallower than `depth`
stay as their own region. Depth is a PRE-SPECIFIED analysis choice: inspect with
`list_depths` (no disease information involved) and fix it before running the signature.
"""
import numpy as np
import pandas as pd


def _paths(ontology):
    out = {}
    for sid, path in zip(ontology["id"], ontology["structure_id_path"]):
        out[int(sid)] = [int(x) for x in str(path).strip("/").split("/") if x]
    return out


def build_region_map(samples, ontology, depth):
    """Return (region_id per sample as array, region table with names)."""
    paths = _paths(ontology)
    unknown = sorted(set(samples["structure_id"].astype(int)) - set(paths))
    if unknown:
        raise ValueError(f"structure_ids not found in ontology: {unknown[:10]}")
    region_ids = np.array([
        (paths[int(s)][depth] if len(paths[int(s)]) > depth else int(s))
        for s in samples["structure_id"]
    ])
    names = ontology.set_index("id")[["acronym", "name"]]
    table = pd.DataFrame({"region_id": sorted(set(region_ids))})
    table["acronym"] = table["region_id"].map(names["acronym"])
    table["name"] = table["region_id"].map(names["name"])
    return region_ids, table


def list_depths(samples_by_donor, ontology, max_depth=10):
    """For each depth: number of regions and sample counts. Contains no disease information."""
    rows = []
    for depth in range(1, max_depth + 1):
        counts = {}
        for samples in samples_by_donor:
            ids, _ = build_region_map(samples, ontology, depth)
            for rid in ids:
                counts[rid] = counts.get(rid, 0) + 1
        values = np.array(list(counts.values()))
        rows.append({"depth": depth, "n_regions": len(counts), "min_samples": values.min(),
                     "median_samples": float(np.median(values)), "max_samples": values.max()})
    return pd.DataFrame(rows)

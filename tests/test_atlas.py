"""All data in this file is SYNTHETIC test data, not Allen Human Brain Atlas content."""
import numpy as np
import pandas as pd

from dbrs.atlas.expression import (differential_stability, regional_gene_matrix, region_means,
                                   select_probe_per_gene, zscore_rows)
from dbrs.atlas.regions import build_region_map


def synthetic_ontology():
    return pd.DataFrame({
        "id": [1, 2, 3, 4, 5],
        "acronym": ["root", "ctx", "ctxA", "ctxB", "sub"],
        "name": ["Root", "Cortex", "Cortex A", "Cortex B", "Subcortex"],
        "structure_id_path": ["/1/", "/1/2/", "/1/2/3/", "/1/2/4/", "/1/5/"],
    })


def test_region_map_depth_and_shallow_structures():
    samples = pd.DataFrame({"structure_id": [3, 4, 5, 2]})
    ids, table = build_region_map(samples, synthetic_ontology(), depth=1)
    assert list(ids) == [2, 2, 5, 2]
    assert set(table["acronym"]) == {"ctx", "sub"}


def test_region_means_and_min_samples():
    expr = pd.DataFrame([[1.0, 3.0, 10.0], [2.0, 4.0, 20.0]], index=["p1", "p2"], columns=[0, 1, 2])
    rm = region_means(expr, [7, 7, 8], ["p1", "p2"], min_samples=2)
    assert list(rm.columns) == [7]
    assert rm.loc["p1", 7] == 2.0


def test_zscore_rows_and_constant_gene_gives_nan():
    df = pd.DataFrame([[1.0, 2.0, 3.0], [5.0, 5.0, 5.0]], index=["g1", "g2"])
    z = zscore_rows(df)
    assert np.allclose(z.loc["g1"], [-1, 0, 1])
    assert z.loc["g2"].isna().all()


def test_probe_selection_prefers_higher_stability():
    genes = pd.Series({"p1": "G", "p2": "G", "p3": "H"})
    ds = pd.Series({"p1": 0.1, "p2": 0.9, "p3": -0.2})
    chosen = select_probe_per_gene(genes, ds)
    assert set(chosen.index) == {"p2", "p3"}


def test_differential_stability_perfect_and_reversed():
    cols = list(range(6))
    a = pd.DataFrame([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]], index=["up", "down"], columns=cols, dtype=float)
    b = pd.DataFrame([[2, 4, 6, 8, 10, 12], [6, 5, 4, 3, 2, 1]], index=["up", "down"], columns=cols, dtype=float)
    ds = differential_stability([a, b], min_shared=5)
    assert np.isclose(ds["up"], 1.0) and np.isclose(ds["down"], -1.0)


def test_regional_gene_matrix_requires_min_donors_and_averages():
    d1 = pd.DataFrame([[1.0, 2.0, 3.0]], index=["G"], columns=[10, 20, 30])
    d2 = pd.DataFrame([[1.0, 2.0, 3.0]], index=["G"], columns=[10, 20, 40])
    z, level, counts = regional_gene_matrix([d1, d2], min_donors=2)
    assert list(z.index) == [10, 20]            # regions 30 and 40 are in only one donor
    assert counts.loc[10] == 2

import json
from dataclasses import dataclass

import pandas as pd

from ..config import data_path

REQUIRED_SAMPLE_COLS = {"structure_id", "structure_acronym", "structure_name", "mni_x", "mni_y", "mni_z"}
REQUIRED_PROBE_COLS = {"probe_id", "gene_symbol", "entrez_id"}


@dataclass
class Donor:
    donor_id: str
    expr: pd.DataFrame      # probes x samples (columns are 0..n-1, aligned with `samples` rows)
    pacall: pd.DataFrame
    samples: pd.DataFrame
    probes: pd.DataFrame
    ontology: pd.DataFrame


def _path(entry, key):
    value = entry[key]
    if isinstance(value, (list, tuple)):
        value = value[0]
    return value


def load_donor(donor_id, entry):
    expr = pd.read_csv(_path(entry, "microarray"), header=None, index_col=0)
    pacall = pd.read_csv(_path(entry, "pacall"), header=None, index_col=0)
    samples = pd.read_csv(_path(entry, "annotation"))
    probes = pd.read_csv(_path(entry, "probes"))
    ontology = pd.read_csv(_path(entry, "ontology"))

    missing = REQUIRED_SAMPLE_COLS - set(samples.columns)
    if missing:
        raise ValueError(f"donor {donor_id}: SampleAnnot is missing columns {sorted(missing)}; has {list(samples.columns)}")
    missing = REQUIRED_PROBE_COLS - set(probes.columns)
    if missing:
        raise ValueError(f"donor {donor_id}: Probes is missing columns {sorted(missing)}; has {list(probes.columns)}")
    if expr.shape[1] != len(samples):
        raise ValueError(f"donor {donor_id}: {expr.shape[1]} expression columns vs {len(samples)} annotated samples")
    if not expr.index.equals(pacall.index):
        raise ValueError(f"donor {donor_id}: expression and PA-call probe order differ")

    expr.columns = range(expr.shape[1])
    pacall.columns = range(pacall.shape[1])
    probes = probes.set_index("probe_id").loc[expr.index]
    samples = samples.reset_index(drop=True)
    return Donor(str(donor_id), expr, pacall, samples, probes, ontology)


def load_all_donors(cfg):
    manifest = json.loads(data_path(cfg, "data", "metadata", "ahba_files.json").read_text())
    wanted = [str(d) for d in cfg["atlas"]["donors"]]
    return [load_donor(d, manifest[d]) for d in wanted]

import json
import platform
from datetime import datetime, timezone
from importlib import metadata

from .config import data_path

TRACKED = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "requests", "abagen"]


def software_versions():
    versions = {"python": platform.python_version()}
    for name in TRACKED:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def record(cfg, step, **details):
    entry = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "step": step,
        **details,
    }
    path = data_path(cfg, "data", "metadata", "provenance.jsonl")
    with open(path, "a") as fh:
        fh.write(json.dumps(entry, default=str) + "\n")
    return entry

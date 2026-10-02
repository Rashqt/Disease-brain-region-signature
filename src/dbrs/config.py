import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_config(path=None):
    path = Path(path) if path else ROOT / "config.toml"
    with open(path, "rb") as fh:
        cfg = tomllib.load(fh)
    cfg["_root"] = ROOT
    return cfg


def data_path(cfg, *parts):
    path = cfg["_root"].joinpath(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path

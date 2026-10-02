"""AHBA microarray download. Uses abagen only as a downloader; analysis code is ours."""
import json

from ..config import data_path


def fetch_ahba(cfg):
    import abagen

    target = data_path(cfg, "data", "raw", "ahba", ".keep").parent
    kwargs = dict(data_dir=str(target), donors=list(cfg["atlas"]["donors"]), verbose=1)
    try:
        files = abagen.fetch_microarray(convert=False, **kwargs)
    except TypeError:  # older abagen without `convert`
        files = abagen.fetch_microarray(**kwargs)
    manifest = {donor: {k: ([str(p) for p in v] if isinstance(v, (list, tuple)) else str(v))
                        for k, v in entry.items()} for donor, entry in files.items()}
    out = data_path(cfg, "data", "metadata", "ahba_files.json")
    out.write_text(json.dumps(manifest, indent=2))
    return manifest

import argparse
import json

import pandas as pd

from .config import load_config


def main(argv=None):
    ap = argparse.ArgumentParser(prog="dbrs", description="Disease Brain-Region Signature")
    ap.add_argument("command", choices=["gwas-inspect", "gwas-fetch", "atlas-fetch", "list-regions",
                                        "atlas-build", "signature", "compare", "all"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--refresh", action="store_true", help="ignore the GWAS cache and re-download")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.command == "gwas-inspect":
        from .gwas.client import extract_records, inspect_first_page
        spec, payload, path = inspect_first_page(cfg)
        recs = extract_records(payload)
        print(f"saved raw first page -> {path}")
        print(f"top-level keys: {list(payload) if isinstance(payload, dict) else 'list'}; records on page: {len(recs)}")
        print("first record:\n" + json.dumps(recs[0], indent=2)[:3000] if recs else "no records")

    elif args.command == "gwas-fetch":
        from .analysis.pipeline import gwas_stage
        assoc, table = gwas_stage(cfg, refresh=args.refresh)
        print(f"{len(assoc)} association-gene rows; {len(table)} genes at p <= {cfg['gwas']['p_threshold']}")

    elif args.command == "atlas-fetch":
        from .atlas.download import fetch_ahba
        fetch_ahba(cfg)
        print("AHBA downloaded; manifest in data/metadata/ahba_files.json")

    elif args.command == "list-regions":
        from .atlas.load import load_all_donors
        from .atlas.regions import list_depths
        donors = load_all_donors(cfg)
        print(list_depths([d.samples for d in donors], donors[0].ontology).to_string(index=False))
        print(f"\ncurrent config region_depth = {cfg['atlas']['region_depth']} (fix this before running the signature)")

    elif args.command == "atlas-build":
        from .analysis.pipeline import atlas_stage
        z, level, meta = atlas_stage(cfg)
        print(f"regional matrix: {z.shape[0]} regions x {z.shape[1]} genes")

    elif args.command in ("signature", "all"):
        from .analysis.pipeline import atlas_stage, gwas_stage, signature_stage
        from .visualization.plots import make_all
        z = level = None
        if args.command == "all":
            gwas_stage(cfg, refresh=args.refresh)
            z, level, _ = atlas_stage(cfg)
        results, null, z, usable = signature_stage(cfg, z, level)
        make_all(cfg, results, null, z, usable)
        cols = ["acronym", "n_samples", "observed", "null_mean", "null_sd", "effect_z", "p_empirical", "p_adj_bh"]
        print(results[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    elif args.command == "compare":
        from .analysis.compare import compare
        results = pd.read_csv(cfg["_root"] / "results" / "signature_results.csv")
        out, source = compare(cfg, results)
        print(f"source: {source}\n{out.to_string(index=False)}")


if __name__ == "__main__":
    main()

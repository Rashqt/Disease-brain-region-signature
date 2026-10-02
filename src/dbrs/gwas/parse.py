"""Turn raw Catalog association records into a flat table, one row per association x gene.

FIELD_CANDIDATES lists the key names we are willing to accept for each logical field.
These are UNVERIFIED guesses for API v2. If a required field cannot be resolved,
SchemaError reports the keys that were actually seen so the mapping can be fixed on purpose.
"""
import re

import pandas as pd

FIELD_CANDIDATES = {
    "association_id": ["association_id", "associationId", "id"],
    "rs_id": ["rs_id", "rsId", "variant_id", "snp_id", "snps"],
    "p_value": ["p_value", "pvalue", "pValue"],
    "study_accession": ["accession_id", "accessionId", "study_accession", "study_id"],
    "pubmed_id": ["pubmed_id", "pubmedId", "pmid"],
    "mapped_genes": ["mapped_genes", "mapped_gene", "mappedGenes"],
    "reported_genes": ["reported_genes", "reported_gene", "author_reported_genes"],
    "disease_trait": ["reported_trait", "disease_trait", "diseaseTrait"],
    "efo_traits": ["efo_traits", "efo_trait", "mapped_trait"],
}
REQUIRED = ["rs_id", "p_value", "mapped_genes"]
NAME_KEYS = ("gene_name", "gene_symbol", "gene", "efo_trait", "trait", "name")
GENE_SPLIT = re.compile(r"\s*[,;]\s*|\s+-\s+")
NOT_A_GENE = {"", "NR", "NA", "NAN", "NONE", "-"}


class SchemaError(RuntimeError):
    pass


def resolve_fields(records):
    seen = set()
    for rec in records:
        seen.update(rec.keys())
    resolved = {}
    for logical, candidates in FIELD_CANDIDATES.items():
        hit = next((c for c in candidates if c in seen), None)
        if hit:
            resolved[logical] = hit
    if "p_value" not in resolved and {"pvalueMantissa", "pvalueExponent"} <= seen:
        resolved["p_value"] = None  # combined from mantissa/exponent
    missing = [f for f in REQUIRED if f not in resolved]
    if missing:
        raise SchemaError(f"Cannot resolve required field(s) {missing}. Keys seen in records: {sorted(seen)}. "
                          "Update FIELD_CANDIDATES in gwas/parse.py.")
    return resolved


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        items = []
        for v in value:
            items.extend(_as_list(v))
        return items
    if isinstance(value, dict):
        for key in NAME_KEYS:
            if isinstance(value.get(key), str):
                return [value[key]]
        return []
    return [str(value)]


def split_genes(value):
    genes = []
    for item in _as_list(value):
        genes.extend(GENE_SPLIT.split(item))
    cleaned = [g.strip().upper() for g in genes]
    return sorted({g for g in cleaned if g not in NOT_A_GENE})


def _p_value(rec, field):
    if field is not None:
        try:
            return float(rec.get(field))
        except (TypeError, ValueError):
            return float("nan")
    try:
        return float(rec["pvalueMantissa"]) * 10.0 ** float(rec["pvalueExponent"])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def parse_associations(records):
    fields = resolve_fields(records)
    get = lambda rec, name: rec.get(fields[name]) if name in fields else None
    rows = []
    for rec in records:
        base = {
            "association_id": get(rec, "association_id"),
            "study_accession": get(rec, "study_accession"),
            "pubmed_id": get(rec, "pubmed_id"),
            "rs_id": ";".join(split_genes_free(get(rec, "rs_id"))),
            "p_value": _p_value(rec, fields["p_value"]),
            "disease_trait": ";".join(_as_list(get(rec, "disease_trait"))),
            "efo_traits": ";".join(_as_list(get(rec, "efo_traits"))),
        }
        for source, name in (("mapped", "mapped_genes"), ("reported", "reported_genes")):
            for gene in split_genes(get(rec, name)):
                rows.append({**base, "gene": gene, "gene_source": source})
    cols = ["association_id", "study_accession", "pubmed_id", "rs_id", "p_value",
            "disease_trait", "efo_traits", "gene", "gene_source"]
    return pd.DataFrame(rows, columns=cols), fields


def split_genes_free(value):
    return [v.strip() for v in _as_list(value) if v.strip()]

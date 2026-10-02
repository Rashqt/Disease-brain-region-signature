import pandas as pd


def build_gene_table(assoc, p_threshold, gene_sources=("mapped",)):
    """Collapse the association x gene table to one row per gene, keeping provenance.

    This does not decide that a gene is causal; it records which associations pointed at it.
    """
    sel = assoc[(assoc["p_value"] <= p_threshold) & assoc["gene_source"].isin(gene_sources)]
    if sel.empty:
        return pd.DataFrame(columns=["gene", "n_associations", "n_variants", "n_studies",
                                     "min_p_value", "gene_sources", "rs_ids", "study_accessions", "pubmed_ids"])
    join = lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x) and str(x) != ""}))
    table = (
        sel.groupby("gene")
        .agg(n_associations=("association_id", "nunique"),
             n_variants=("rs_id", "nunique"),
             n_studies=("study_accession", "nunique"),
             min_p_value=("p_value", "min"),
             gene_sources=("gene_source", join),
             rs_ids=("rs_id", join),
             study_accessions=("study_accession", join),
             pubmed_ids=("pubmed_id", join))
        .reset_index()
        .sort_values("gene")
        .reset_index(drop=True)
    )
    return table

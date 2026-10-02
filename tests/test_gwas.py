"""All data in this file is SYNTHETIC test data, not real GWAS Catalog content."""
import pandas as pd
import pytest

from dbrs.gwas.genes import build_gene_table
from dbrs.gwas.parse import SchemaError, parse_associations, split_genes


def synthetic_records():
    return [
        {"association_id": 1, "rs_id": "rsTEST1", "p_value": 1e-10, "accession_id": "GCSTX1", "pubmed_id": "1",
         "mapped_genes": ["GENEA", "GENEB"], "reported_genes": "NR"},
        {"association_id": 2, "rs_id": "rsTEST2", "p_value": 1e-9, "accession_id": "GCSTX2", "pubmed_id": "2",
         "mapped_genes": "GENEA, GENEC", "reported_genes": ["GENED"]},
        {"association_id": 3, "rs_id": "rsTEST3", "p_value": 1e-3, "accession_id": "GCSTX2", "pubmed_id": "2",
         "mapped_genes": ["GENEZ"], "reported_genes": []},
        {"association_id": 4, "rs_id": "rsTEST4", "p_value": 2e-12, "accession_id": "GCSTX3", "pubmed_id": "3",
         "mapped_genes": "HLA-DRB1 - GENEE", "reported_genes": None},
    ]


def test_split_genes_handles_lists_strings_and_placeholders():
    assert split_genes("a, b;NR") == ["A", "B"]
    assert split_genes(["x", None, "x"]) == ["X"]
    assert split_genes("HLA-DRB1 - GENEE") == ["GENEE", "HLA-DRB1"]  # hyphen inside a symbol is not a separator


def test_parse_keeps_mapped_and_reported_separate():
    df, _ = parse_associations(synthetic_records())
    assert set(df.loc[df.gene == "GENED", "gene_source"]) == {"reported"}
    assert set(df.loc[df.gene == "GENEA", "gene_source"]) == {"mapped"}


def test_gene_table_threshold_duplicates_and_provenance():
    df, _ = parse_associations(synthetic_records())
    table = build_gene_table(df, 5e-8, ("mapped",)).set_index("gene")
    assert "GENEZ" not in table.index            # p too weak
    assert "GENED" not in table.index            # reported-only genes excluded by default
    assert table.loc["GENEA", "n_associations"] == 2
    assert table.loc["GENEA", "study_accessions"] == "GCSTX1;GCSTX2"


def test_missing_required_field_fails_loudly():
    with pytest.raises(SchemaError):
        parse_associations([{"foo": 1}])


def test_mantissa_exponent_p_value():
    recs = [{"rs_id": "rsT", "pvalueMantissa": 2, "pvalueExponent": -9, "mapped_genes": ["G1"]}]
    df, _ = parse_associations(recs)
    assert df.loc[0, "p_value"] == pytest.approx(2e-9)

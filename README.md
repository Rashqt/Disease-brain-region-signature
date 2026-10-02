# Disease Brain-Region Signature

Do the genes that GWAS link to a neurological disease show a distinctive expression pattern across human brain regions, and does that pattern line up with the regions the disease is known to hit?

That's the whole question. This is an investigation, not a claim that I'm finding the cause of anything. A boring null result is a fine outcome and I'll report it as one.

> **Status: Stage 1 (v0.1).** The code is written and the unit tests pass on synthetic data, but the pipeline has **not been run on the real GWAS Catalog / Allen data yet**. No results exist. Nothing here should be read as a finding.

## Why I built this

I code, I like biology, and I wanted a project where the stats have to be right or the whole thing is meaningless. GWAS hits are easy to over-read and brain heatmaps are easy to make look convincing, so the point is to build the careful version: matched null models, multiple-testing correction, and a hard wall between the computed result and what's already known about the disease.

## First test case: Parkinson's disease

Chosen for data suitability, not for how nice the output might look:

- large, well-replicated GWAS with many independent loci, so the gene set is big enough for a null model to mean something
- the Allen Human Brain Atlas (AHBA) samples subcortical structures, not just cortex
- Alzheimer's signal is dominated by microglia/immune biology that bulk regional expression handles poorly; Huntington's is monogenic; ALS vulnerability sits partly in spinal cord, which AHBA doesn't cover

Known weak spot: substantia nigra is sparsely sampled in AHBA, so the region most associated with PD may be underpowered. That will be documented, not hidden.

## Pipeline

```
GWAS Catalog REST API v2  ->  associations (p <= 5e-8)  ->  Catalog-mapped genes (provenance kept)
AHBA microarray (6 donors) ->  probe filter -> best probe per gene -> region x gene z-scores
disease genes + atlas      ->  observed regional score vs expression-matched random gene sets
                           ->  empirical p, effect size, bootstrap CI, BH-FDR across regions
                           ->  figures
(only afterwards)          ->  compare with independently cited vulnerable regions
```

### Methods in plain words

- **Genes from GWAS.** Variants are not genes. I use the genes the Catalog itself maps to each variant (`mapped`), keep author-reported genes in the tables but out of the default set, and store rsID, study accession and PMID for every gene. A gene in the set means "a genome-wide-significant association points near it", not "it causes the disease".
- **Gene matching to the atlas.** Exact symbol match. Genes that don't match are listed in `gene_set_mapping_report.csv`, not silently dropped or remapped.
- **Regions.** Allen structure ontology, collapsed to a fixed depth (`region_depth` in `config.toml`). `list-regions` shows the options and contains no disease info, so the depth gets fixed before any result is seen.
- **Expression.** Probes need a present call in >= 50% of samples; per gene the probe with the highest differential stability across donors is kept. Per donor, each gene is z-scored across regions (so heavily sampled cortex doesn't dominate), then averaged over donors. Regions need >= 3 donors.
- **Score and null.** Region score = mean z of the disease genes. Null = 10,000 random gene sets of the same size, each gene matched on expression-level bin, excluding disease genes. Empirical one-sided p (enrichment, pre-specified), add-one so it's never 0, BH-FDR across regions, effect as (obs - null mean)/null SD, plus a gene-bootstrap 95% CI.
- **Known vulnerable regions.** `data/metadata/known_vulnerable_regions.toml`, filled from a cited source before `compare` will run. Never used to tune anything.

## Install and run

Python 3.11+.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
pytest
```

```bash
python -m dbrs.cli gwas-inspect    # look at one raw API page first (the v2 schema is not verified yet)
python -m dbrs.cli gwas-fetch
python -m dbrs.cli atlas-fetch     # ~4 GB download
python -m dbrs.cli list-regions    # choose region_depth in config.toml, then stop touching it
python -m dbrs.cli atlas-build
python -m dbrs.cli signature
python -m dbrs.cli compare         # needs the cited known-regions file
```

AHBA data is non-commercial-use and not redistributed here. `atlas-fetch` downloads it from the Allen Institute.

## Layout

```
config.toml            every analysis parameter + seed
src/dbrs/gwas/         API client, parsing, gene table
src/dbrs/atlas/        download, loading, region mapping, expression
src/dbrs/stats/        matched null, p-values, BH
src/dbrs/analysis/     pipeline stages, known-region comparison
src/dbrs/visualization/
data/{raw,processed,metadata}/   raw/processed are git-ignored; metadata/provenance.jsonl logs every run
results/               tables + figures
tests/                 synthetic data only, labelled as such
```

## What works / what's experimental

- Works (unit-tested on synthetic data): gene parsing and dedup, region mapping, probe selection, aggregation, null sampling, p-values, BH, seed reproducibility.
- Untested on real data: everything that touches the network or the real files.
- Experimental: the GWAS v2 field names and pagination in `gwas/parse.py` and `config.toml` are best guesses. The code fails loudly if they're wrong rather than guessing.

## Limitations (so far)

- Catalog gene mapping is positional (overlap / nearest), not functional. Many mapped genes won't be the causal gene at their locus.
- Nearby genes at one locus have correlated expression; the null draws independent genes, which can make p-values optimistic.
- Null matches expression level only. Gene length, GC content and others are not matched yet.
- Healthy adult bulk tissue, 6 donors (5 male), microarray. No cell types, no disease tissue, no age effects.
- Allen MNI coordinates have known inaccuracies; the scatter plot is a sample-location view, not a surface map.
- The PD trait in the Catalog may include sub-phenotypes (e.g. onset age); trait filtering is a Stage 5 decision.
- Hemispheres are pooled.

## Roadmap

Run it and fix what breaks (2-3), trait filtering + sensitivity analyses + extra confounder matching (5), then maybe more diseases. ML only if it earns its place.

Coded with late night coffee, questionable sleep, and too many datasets. - Rashqt

MIT license.

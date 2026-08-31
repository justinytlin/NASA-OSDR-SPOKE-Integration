# NASA-OSDR-SPOKE-Integration

Pipelines for projecting NASA OSDR / GeneLab omics datasets onto the
[SPOKE](https://spoke.ucsf.edu/) biomedical knowledge graph.

Two complementary routes, from a differential-expression table to ranked
SPOKE nodes (diseases, pathways, compounds, symptoms, ...):

1. **Full-graph PSEV embeddings** (`replication/`) — the method of
   [Nelson et al., *Life* 2021, 11:42](https://doi.org/10.3390/life11010042),
   replicated and generalized: weight precomputed per-gene PSEVs
   (Propagated SPOKE Entry Vectors) by log2 fold change and rank all
   389,297 SPOKE v2 nodes. Deep multi-hop signal; needs a one-time 25.6 GB
   download.
2. **Live REST-API projections** (`spoke-api-client/`) — lightweight client
   for the public SPOKE REST API: neighborhood projections, subgraph
   fetching, and PSEV-style personalized-PageRank per-sample embeddings.
   No bulk data needed; current graph instead of the frozen 2019 one.

## Layout

```
replication/
  psev_pipeline.py     GeneLab FC table -> PSEV node ranks (faithful port of
                       SPOKE_NASA genelab_fc_to_psev.py; exposes
                       compute_rank_table() for reuse)
  osdr_to_psev.py      ANY OSDR/GeneLab DE table (CSV or TSV, mouse or
                       --human) -> full-SPOKE rank table
  run_all.py           batch run of the paper's six GLDS studies
  welch_meta.py        multi-study Welch's t-test meta-analysis (Methods 2.5)
  plot_violins.py      violin plots of rank distributions
  test_synthetic.py    end-to-end smoke test against a small fake PSEV zip
spoke-api-client/
  spoke.py             CLI for the public SPOKE REST API (see its README)
  spoke_api/           client library + MGI/HomoloGene ortholog mapping
  examples/
    osdr_spoke_projection.py   DE table -> 1-hop disease/pathway hit counts
    fetch_subgraph.py          union of 1-hop neighborhoods -> subgraph.jsonl
    spoke_embeddings_ml.py     per-sample PPR embeddings + LOO classification
    welch_nodes.py             which nodes separate two sample groups
                               (Welch per node, rank-by-type, top-2.5%)
    predict_treatments.py, alzheimers_multidb.py, disease_gene_table.py
psev_data/             (gitignored) PSEV matrices + ortholog tables, see below
SPOKE_NASA/            (not tracked) clone of baranzini-lab/SPOKE_NASA
```

## Data prerequisites

Only needed for the full-graph PSEV route (`replication/`):

| What | Where | Put it at |
|---|---|---|
| SPOKE v2 node universe + paper ground truth | `git clone https://github.com/baranzini-lab/SPOKE_NASA` | `SPOKE_NASA/` |
| Precomputed gene PSEVs (25.6 GB) | Zenodo [10.5281/zenodo.4404618](https://doi.org/10.5281/zenodo.4404618) (`psev_data/parallel_download.sh` helps) | `psev_data/gene_psev.zip` |
| Frozen NCBI HomoloGene build 68 | `https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data` | `psev_data/homologene_build68.data` |

For the API route, optionally the MGI ortholog report
(`HOM_MouseHumanSequence.rpt` from https://www.informatics.jax.org/downloads/reports/)
for `--orthologs`. Note: modern MGI reports dropped the HomoloGene ID column,
which is why the PSEV route pins build 68.

Python 3.10+, `numpy pandas scipy requests`; `scikit-learn matplotlib openpyxl`
for the ML/meta extras.

## Quickstart: any OSDR study -> SPOKE node ranks

```bash
# grab a processed DE table (CSV) from OSDR, e.g. via
#   https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-NNN/files/
#   https://osdr.nasa.gov/geode-py/ws/studies/OSD-NNN/download?source=datamanager&file=...

# check gene mapping + contrast selection (free, no big data touched)
python3 replication/osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01 --prep-only

# full run (~10-15 min: two streaming passes over gene_psev.zip)
python3 replication/osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01
```

Output: `replication/psev_out/<name>_ranks_and_rank_by_type_for_meta0_1.tsv` —
one row per SPOKE node with per-contrast overall rank and `Rank_by_type_*`
columns (rank 1 = strongest within that node type). The format matches the
Nelson-2021 study outputs, so `welch_meta.py`-style pooling works across old
and new studies. Human datasets: add `--human`. Gene filters:
`--filter none | pvalue | same-direction` (the last is the paper's rule for
space/ground/basal designs).

## Quickstart: which SPOKE nodes separate two sample groups

```bash
python3 spoke-api-client/examples/fetch_subgraph.py ...      # signature genes -> subgraph.jsonl
python3 spoke-api-client/examples/spoke_embeddings_ml.py ... # per-sample PPR embeddings
python3 spoke-api-client/examples/welch_nodes.py --dir OUTDIR
```

Treat the resulting per-node p-values as rankings, not calibrated
significance (embedding dimensions are correlated and n is small).

## Testing

`python3 replication/test_synthetic.py [SCRATCH_DIR]` builds a small fake
PSEV zip and runs the pipeline + meta-analysis end to end (requires
`SPOKE_NASA/` and the HomoloGene file, not the 25.6 GB zip).

# NASA-OSDR-SPOKE-Integration

1. **[The PSEV pipeline](#1-charlottes-paper-replication-the-psev-pipeline)** —
   a replication of Charlotte Nelson's spaceflown-mice paper, generalized so
   any OSDR dataset can be run through it. Needs a one-time 25.6 GB download and uses
   a frozen 2019 snapshot of SPOKE.
3. **[The public SPOKE API](#2-the-public-spoke-api-method)** — a lightweight
   client for the live SPOKE REST API. No bulk data, always-current graph,
   but limited to local neighborhoods around your genes.
4. **[Anatomy sanity check](#3-sanity-check-do-the-anatomy-nodes-make-sense)** —
   every dataset comes from a known tissue, so the tissue is free ground
   truth: if the rankings are working, the source tissue should appear ranked highly
   extreme of the Anatomy ranking.

## The overall workflow for a new dataset

```
1. Download the differential-expression CSV for OSD-NNN from OSDR
2. replication/osdr_to_psev.py --prep-only     check gene mapping (free)
3. replication/osdr_to_psev.py                 full run -> ranks for all 389,297 nodes
4. replication/osdr_null_control.py            ground-control-only null -> which
                                               top nodes beat a no-biology control
5. replication/anatomy_check.py                does the source tissue surface?
6. Interpret: nodes in the top 2.5% of their type that beat all null splits
   (welch_meta.py can pool several studies into a meta-analysis)
```

Requirements: Python 3.10+ with `numpy pandas scipy requests`
(`scikit-learn matplotlib openpyxl` for the ML and meta-analysis extras).

---

## 1. Charlotte's paper replication (the PSEV pipeline)

[Nelson et al., *Life* 2021, 11:42](https://doi.org/10.3390/life11010042)
embedded six GeneLab mouse transcriptomes (thymus, spleen, liver; GLDS-4,
-244, -245, -246, -288, -289) into SPOKE and found spaceflight-associated
nodes. The method: every gene in SPOKE has a precomputed **PSEV** (Propagated
SPOKE Entry Vector) — a random-walk profile of how strongly that gene reaches
every node in the graph. Weight the signature genes' PSEVs by their log2 fold
changes, sum, and rank all 389,297 nodes; then rank *within each node type*
so "top disease" isn't drowned out by the 287k compounds.

### Files

| File | What it does |
|---|---|
| `replication/psev_pipeline.py` | The core algorithm (faithful port of the paper's `genelab_fc_to_psev.py`); streams the 25.6 GB PSEV zip in two passes so it never needs 60 GB of disk. Exposes `compute_rank_table()` for reuse. |
| `replication/run_all.py` | Batch-runs the paper's six studies. |
| `replication/welch_meta.py` | The paper's meta-analysis: pools comparisons into Space-v-Ground / Space-v-Basal / Ground-v-Basal groups, Welch's t-test per node, top 2.5% per type, and self-scores against the paper's published node list. |
| `replication/osdr_to_psev.py` | **The generalization** — runs ANY OSDR/GeneLab DE table (CSV or TSV, mouse or `--human`) through the same machinery, producing output in the same format so studies old and new can be pooled. |
| `replication/osdr_null_control.py` | **The control** — splits the ground-control samples into halves (no biological difference), runs the fake contrasts through the identical machinery, and flags which real findings fail to beat this no-biology null (hubs and baseline tissue variability rank high for *any* gene set). |
| `replication/plot_violins.py` | Rank-distribution figures. |
| `replication/test_synthetic.py` | End-to-end smoke test against a small fake PSEV zip (no 25.6 GB download needed). |

### One-time data setup

| What | Where | Put it at |
|---|---|---|
| SPOKE v2 node universe + paper ground truth | `git clone https://github.com/baranzini-lab/SPOKE_NASA` | `SPOKE_NASA/` |
| Precomputed gene PSEVs (25.6 GB) | Zenodo [10.5281/zenodo.4404618](https://doi.org/10.5281/zenodo.4404618) (`psev_data/parallel_download.sh` helps) | `psev_data/gene_psev.zip` |
| Frozen NCBI HomoloGene build 68 | `https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data` | `psev_data/homologene_build68.data` |

(Build 68 is pinned because modern MGI ortholog reports dropped the
HomoloGene ID column the paper relied on.)

### Running a new dataset

```bash
# DE tables come from OSDR, e.g. list files at
#   https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-NNN/files/
# then download via
#   https://osdr.nasa.gov/geode-py/ws/studies/OSD-NNN/download?source=datamanager&file=...

# 1. check gene mapping and contrast selection (free)
python3 replication/osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01 --prep-only

# 2. full run (~10-15 min: two streaming passes over the zip)
python3 replication/osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01

# 3. ground-control null (all splits share one pass over the zip)
python3 replication/osdr_null_control.py TABLE.csv --name OSD-NNN-GCnull \
    --sample-pattern _GC_ --n-genes <signature size from step 2> \
    --real replication/psev_out/OSD-NNN_ranks_and_rank_by_type_for_meta0_1.tsv \
    --real-col "Log2fc_(Space Flight)v(Ground Control)"
```

### Reading the output

`psev_out/<name>_ranks_and_rank_by_type_for_meta0_1.tsv` has one row per
SPOKE node. Use the `Rank_by_type_*` columns: **rank 1 = most associated
with UP-regulated genes** of that node type, rank ≈ type size = most
associated with DOWN-regulated genes. Both tails are findings; the middle is
noise. "High" conventionally means the top 2.5% of the node's type.

Single-study ranks carry no significance — that's what the null control
adds. Its `<name>_vs_real.tsv` gives each node an empirical percentile
against the no-biology splits (`p_top`, floor 1/(n_splits+1)) and a
`beats_all_nulls` flag. Report nodes that are top-2.5% **and** beat all
nulls; a node the null also produces (e.g. "carcinoma") is a
graph-connectivity artifact, not biology. Gene filters: `--filter none |
pvalue | same-direction` (the last is the paper's rule for
space/ground/basal designs). Human datasets: add `--human`.

---

## 2. The public SPOKE API method

The [public SPOKE REST API](https://spoke.rbvi.ucsf.edu/swagger/) needs no
authentication and reflects the *current* graph, but only serves
neighborhoods — so this route builds small subgraphs around your genes
instead of using the whole graph. Useful as an independent check on the PSEV
route and for anything that needs up-to-date graph content.

| File | What it does |
|---|---|
| `spoke-api-client/spoke.py` + `spoke_api/` | CLI and library: search, node lookup, neighborhoods, SEA compound similarity; ortholog mapping helpers. See `spoke-api-client/README.md` for every command. |
| `examples/osdr_spoke_projection.py` | Simplest analysis: signature genes → 1-hop neighborhoods → which diseases/pathways get hit. Fast but hub-biased. |
| `examples/fetch_subgraph.py` | Unions the 1-hop neighborhoods of the signature genes into a local subgraph (`subgraph.jsonl`). |
| `examples/spoke_embeddings_ml.py` | PSEV-style per-sample embeddings on that subgraph (personalized PageRank on up/down-regulated genes) + leave-one-out classification. |
| `examples/welch_nodes.py` | Which nodes *separate* two sample groups: Welch's t-test per node on the embeddings, rank-by-type, top-2.5% flags. Turns "the classes separate" into "these nodes separate them". |
| `examples/predict_treatments.py`, `alzheimers_multidb.py`, `disease_gene_table.py` | Drug-repurposing-style queries on the live graph. |

```bash
python3 spoke-api-client/examples/fetch_subgraph.py ...      # genes -> subgraph.jsonl
python3 spoke-api-client/examples/spoke_embeddings_ml.py ... # per-sample embeddings
python3 spoke-api-client/examples/welch_nodes.py --dir OUTDIR
```

Treat `welch_nodes.py` p-values as rankings, not calibrated significance
(embedding dimensions are correlated and sample counts are small). When
fetching subgraphs, include `Anatomy` in the node filters if you want the
anatomy sanity check to work on this route.

---

## 3. Sanity check: do the anatomy nodes make sense?

Every dataset's source tissue is known ground truth, which makes the Anatomy
ranking a free positive control for the whole pipeline:

```bash
python3 replication/anatomy_check.py RANKS.tsv --tissue "spleen"
```

It prints the top and bottom Anatomy nodes and the source tissue's
percentile. **Interpret direction, not just position**: the ranking reflects
the *differential* signature, so a tissue whose identity genes are
suppressed by the condition correctly lands at the *bottom* — e.g. under
spaceflight immune suppression, spleen ranks in the bottom 2% of its own
spleen dataset, while a hippocampus dataset puts hippocampal/basal-forebrain
structures near the top. Either extreme is a pass; the source tissue sitting
mid-ranking is the red flag.

By default the check hides UBERON life-stage/fluid pseudo-anatomy nodes
("pupal stage", "colostrum") whose sparse gene annotations make them
unstable at the extremes (`--include-stages` to show them), and it pools
Space-v-Ground-style comparisons automatically (`--col` / `--col-regex` for
other designs). Also watch for near-duplicate node families (ten thoracic
dorsal-root-ganglion nodes = one finding).

---

## Testing

`python3 replication/test_synthetic.py [SCRATCH_DIR]` fabricates a small
random PSEV zip over the real node universe and runs the pipeline plus
meta-analysis end to end in minutes — the answer to "how do you test a
pipeline whose input is 25.6 GB".

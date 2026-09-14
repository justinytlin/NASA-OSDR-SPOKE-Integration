# replication/ — the PSEV pipeline

The replication of [Nelson et al., *Life* 2021, 11:42](https://doi.org/10.3390/life11010042)
(six GeneLab mouse transcriptomes embedded into SPOKE via Propagated SPOKE
Entry Vectors), generalized so any OSDR dataset can be run through it.

The method in one paragraph: every gene in SPOKE has a precomputed **PSEV**,
a random-walk profile of how strongly that gene reaches every node in the
graph. Weight the signature genes' PSEVs by their log2 fold changes, sum, and
rank all nodes; then rank *within each node type* so "top disease" isn't
drowned out by hundreds of thousands of compounds.

## Scripts

| Script | Role | What it does |
|---|---|---|
| `psev_pipeline.py` | core | Faithful port of the paper's `genelab_fc_to_psev.py`. Streams the 25.6 GB PSEV zip (or a `psev_build/` directory) in two passes so it never needs 60 GB of disk. Exposes `compute_rank_table()` and `PsevStore`. |
| `osdr_to_psev.py` | **DE route** | Any OSDR/GeneLab differential-expression table (CSV or TSV, mouse or `--human`) → one rank table per contrast, same format as the paper's outputs. |
| `counts_to_psev.py` | **counts route** | Unnormalized counts CSV → one PSEV embedding *per sample*, weighted by log2(CPM+1) or per-gene z-scores. An extension of the paper, not part of the replication. |
| `osdr_null_control.py` | control | Splits ground-control samples into halves (no biology), runs the fake contrasts through the identical machinery, and flags which real findings fail to beat this null. |
| `anatomy_check.py` | sanity check | Top/bottom Anatomy nodes, source-tissue percentile, and `--tissue-regex` family scoring with a Mann-Whitney test. Works on DE-route and per-sample columns. |
| `run_all.py` | replication | Batch-runs the paper's six studies. |
| `welch_meta.py` | replication | The paper's meta-analysis: pools comparisons, Welch's t-test per node, top 2.5% per type, self-scores against the published node list. `--spoke` points it at another node universe. |
| `plot_violins.py`, `plot_anatomy_families.py` | figures | Rank-distribution violins; anatomy-family strip plots across one or more PSEV builds. |
| `test_synthetic.py` | test | Fabricates a small random PSEV zip over the real node universe and runs pipeline + meta-analysis end to end in minutes. |

Requirements: Python 3.10+ with `numpy pandas scipy` (`matplotlib` for the
figures, `openpyxl` for the meta-analysis ground-truth file).

## One-time data setup

| What | Where | Put it at |
|---|---|---|
| SPOKE v2 node universe + paper ground truth | `git clone https://github.com/baranzini-lab/SPOKE_NASA` | `../SPOKE_NASA/` |
| Precomputed gene PSEVs (25.6 GB) | Zenodo [10.5281/zenodo.4404618](https://doi.org/10.5281/zenodo.4404618) (`../psev_data/parallel_download.sh` helps) | `../psev_data/gene_psev.zip` |
| Frozen NCBI HomoloGene build 68 | `https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data` | `../psev_data/homologene_build68.data` |

Build 68 is pinned because modern MGI ortholog reports dropped the HomoloGene
ID column the paper relied on. Paths are relative to the repository root; the
scripts default to these locations and every one takes `--zip`, `--spoke`
and `--homologene` overrides. A PSEV built with `../psev_build/` is used by
passing its directory as both `--zip` and `--spoke`.

## DE route: a new dataset in three commands

```bash
# DE tables come from OSDR, e.g. list files at
#   https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-NNN/files/
# then download via
#   https://osdr.nasa.gov/geode-py/ws/studies/OSD-NNN/download?source=datamanager&file=...

# 1. check gene mapping and contrast selection (free)
python3 osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01 --prep-only

# 2. full run (~10-15 min: two streaming passes over the zip)
python3 osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01

# 3. ground-control null (all splits share one pass over the zip)
python3 osdr_null_control.py TABLE.csv --name OSD-NNN-GCnull \
    --sample-pattern _GC_ --n-genes <signature size from step 2> \
    --real psev_out/OSD-NNN_ranks_and_rank_by_type_for_meta0_1.tsv \
    --real-col "Log2fc_(Space Flight)v(Ground Control)"
```

Gene filters: `--filter none | pvalue | same-direction` (the last is the
paper's rule for space/ground/basal designs). Human datasets: add `--human`.

## Reading a rank table

`psev_out/<name>_ranks_and_rank_by_type_for_meta0_1.tsv` has one row per
SPOKE node and, per contrast (or per sample on the counts route), two
columns:

- the overall rank across all nodes (high = most associated). Dominated by
  node type: genes and proteins sit high because genes are where the walks
  start. Rarely what you want.
- `Rank_by_type_<contrast>`: the same score ranked within the node's type.
  **Rank 1 = most associated with UP-regulated genes**; rank ≈ type size =
  most associated with DOWN-regulated genes. Divide by the type size for a
  percentile. Both tails are findings; the middle is noise. "High"
  conventionally means the top 2.5% of the type.

Single-study ranks carry no significance; that is what the null control
adds. Its `<name>_vs_real.tsv` gives each node an empirical percentile
against the no-biology splits (`p_top`, floor 1/(n_splits+1)) and a
`beats_all_nulls` flag. Report nodes that are top-2.5% **and** beat all
nulls. A node the null also produces (e.g. "carcinoma") is a
graph-connectivity artifact, not biology. Low-degree nodes swing to the
extremes for any gene set; look up a node's degree before believing its rank.

## Counts route: one embedding per sample

```bash
# counts come from the same OSDR file listing, e.g.
#   .../OSD-564/download?source=datamanager&file=GLDS-569_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv

# check the ENSEMBL -> Entrez -> human mapping and the sample selection (free)
python3 counts_to_psev.py COUNTS.csv --name OSD-NNN-GC --samples _GC_ --gene-map DE_TABLE.csv --prep-only

# full run: log2(CPM+1)-weighted PSEV per kept sample
python3 counts_to_psev.py COUNTS.csv --name OSD-NNN-GC --samples _GC_ --gene-map DE_TABLE.csv [--weight zscore] [--min-cpm 1]
```

`--samples` is a substring filter on column names (repeatable); `--gene-map`
is any GeneLab DE table for the same organism, used only for its
ENSEMBL/ENTREZID columns. Expression-weighted profiles of samples from one
tissue are near-identical to each other (rank correlation > 0.999), so use a
contrast or `--weight zscore` when the question is *differences between
samples*. `../gc_counts_check/` is the systematic version of this route
across many OSDR datasets.

## Anatomy sanity check

Every dataset's source tissue is known ground truth, which makes the Anatomy
ranking a free positive control:

```bash
python3 anatomy_check.py RANKS.tsv --tissue "spleen"                       # DE route
python3 anatomy_check.py RANKS.tsv --col-regex "_GC_" \
    --tissue-regex "hippocamp|dentate gyrus|Ammon|subiculum"               # counts route, node family
```

**Interpret direction, not just position**: the DE-route ranking reflects
the *differential* signature, so a tissue whose identity genes are suppressed
by the condition correctly lands at the *bottom*. Either extreme is a pass;
the source tissue sitting mid-ranking is the red flag. The check hides UBERON
life-stage/fluid pseudo-anatomy nodes by default (`--include-stages` shows
them) and pools Space-v-Ground-style comparisons automatically (`--col` /
`--col-regex` for other designs). Near-duplicate node families (ten thoracic
dorsal-root-ganglion nodes) are one finding, not ten.

This check is informative on the 2019 PSEVs and **not** on builds from the
2025 SPOKE snapshot: see `../psev_build/README.md` ("Anatomy-edge variants")
and `figures/osd564_gc_hippocampus_anatomy_four_builds.png`.

## figures/

`osd564_gc_hippocampus_anatomy*.png|svg`: OSD-564 ground-control expression
profiles, where hippocampal nodes sit among Anatomy nodes, on the 2019 PSEVs
alone, 2019 vs 2025, and across four builds. Produced by
`plot_anatomy_families.py`.

## Testing

`python3 test_synthetic.py [SCRATCH_DIR]` builds a small random PSEV zip over
the real node universe and runs the pipeline plus meta-analysis end to end in
minutes, the answer to "how do you test a pipeline whose input is 25.6 GB".

# Worked example: OSD-564 (mouse hippocampus, Rodent Research-10)

A complete run of the PSEV workflow on
[OSD-564](https://osdr.nasa.gov/bio/repo/data/studies/OSD-564): right
hippocampus of mice flown on the ISS (RR-10 / SpaceX-21), bulk RNA-seq,
Space Flight vs Ground Control. This directory holds the real outputs plus
this walkthrough of how to read them.

## How these files were produced

```bash
# DE table: GLDS-569_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv
# from https://osdr.nasa.gov/geode-py/ws/studies/OSD-564/download?...

# 1. main run: 1,415 genes at p<0.01 -> ranks for all 389,297 SPOKE nodes
python3 replication/osdr_to_psev.py OSD-564_table.csv --name OSD-564 \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01

# 2. ground-control null: all 10 possible 2-vs-3 splits of the 5 GC samples
python3 replication/osdr_null_control.py OSD-564_table.csv --name OSD-564-GCnull \
    --sample-pattern _GC_ --n-genes 1415 \
    --real replication/psev_out/OSD-564_ranks_and_rank_by_type_for_meta0_1.tsv \
    --real-col "Log2fc_(Space Flight)v(Ground Control)"

# 3. anatomy sanity check
python3 replication/anatomy_check.py \
    replication/psev_out/OSD-564_ranks_and_rank_by_type_for_meta0_1.tsv \
    --tissue "Ammon's horn"
```

## What's here

| File | What it is |
|---|---|
| `OSD-564_ranks_..._meta0_1.tsv.gz` | Main result: every SPOKE node with its overall rank and `Rank_by_type_` rank for the Space-v-Ground contrast. (`pandas.read_csv(..., sep="\t")` reads .gz directly.) |
| `OSD-564-GCnull_vs_real.tsv.gz` | Null-control comparison: per node, real rank vs the 10 no-biology splits (`null_median`, `null_best`, `p_top`, `beats_all_nulls`). |
| `top_findings.csv` | The distilled answer: top 25 robust nodes per node type (top 2.5% of their type AND beating every null split). |
| `anatomy_check_output.txt` | The tissue sanity check as run. |
| `run_psev.log.txt`, `run_null.log.txt` | Logs of both runs. |

(The full null rank ensemble is ~70 MB and is summarized by the `_vs_real`
file; regenerate it with the command above if you need per-split ranks.)

## How to interpret the results, step by step

**1. Start from `Rank_by_type_Log2fc_(Space Flight)v(Ground Control)`** in
the main table. Rank 1 = the node most associated with genes UP-regulated in
spaceflight, within its node type; rank ≈ type size = most associated with
DOWN-regulated genes. Both tails are findings — the middle is noise. The
top-2.5% convention: `rank / type_size <= 0.025`.

**2. Check the biology passes the smell test** (`anatomy_check_output.txt`).
For a hippocampus dataset, the up tail is basal forebrain circuitry
(diagonal band of Broca, stria terminalis; Ammon's horn sits at the 15.8th
percentile), and the down tail is lymphoid tissue — the well-documented
spaceflight immune suppression showing up even in brain. The source tissue
at an *extreme* (either one) is a pass; source tissue mid-ranking would be
the red flag.

**3. Filter through the null control** (`_vs_real.tsv.gz`). A single-study
rank has no significance; the null asks "does a no-biology contrast built
from ground-control samples alone also produce this node?" Two instructive
rows:

| Node | real rank | null best | verdict |
|---|---|---|---|
| carcinoma (Disease) | 1 | 1 | a null split also ranks it #1 → hub artifact, discard |
| diagonal band of Broca (Anatomy) | 1 | 138 | no null comes close → robust finding |

`p_top` is the empirical percentile; with 10 splits its floor is 1/11 ≈ 0.09,
so read `beats_all_nulls` as a robustness flag, not a p-value. Of the 9,739
top-2.5% nodes, 8,711 (89%) beat all nulls — the other 11% (including the
generic cancer hubs atop the Disease list) are connectivity artifacts.

**4. Read `top_findings.csv`** — that filtering already applied, top 25 per
node type. Highlights that emerge: interneuron-migration / GABAergic
development processes, NMDA-receptor and axon-guidance pathways, basal
forebrain anatomy, and hypothalamic-flavored symptoms (hypothermia,
hyperphagia). Two nodes from the Nelson paper's showcase replicate here:
Vitamin D (calciferol) metabolism (top 2.2% of pathways) and regulation of
cortisol secretion (top 1.4% of biological processes).

**5. Standard caveats.** Ranks come from the frozen 2019 SPOKE v2 graph;
near-duplicate nodes (same pathway from two sources) should be counted once;
top Compounds are mostly obscure ChEMBL screening molecules — filter to
approved drugs before drawing countermeasure conclusions; and one study is
hypothesis generation — pooling several studies (`welch_meta.py`) is what
buys real p-values.

## The same dataset through the public SPOKE API (`api_projection/`)

For comparison, `api_projection/` holds OSD-564 run through the *other*
method — `spoke-api-client/examples/osdr_spoke_projection.py`, which queries
the live SPOKE REST API and counts which diseases/pathways the signature
genes' 1-hop neighborhoods hit:

| File | What it is |
|---|---|
| `gene_matches.csv` | The ~100 top signature genes, their human orthologs, and whether SPOKE knows them. |
| `diseases.csv`, `pathways.csv` | Diseases/pathways ranked by how many signature genes touch them directly. |
| `edges.csv` | Every gene→disease / gene→pathway edge found (for graph viz). |

Putting the two methods side by side is instructive. The API route's top
diseases are **generic hubs** — "central nervous system disease", "disease",
"cancer" — because raw hit-counting rewards nodes with many gene edges and
sees only one hop. The PSEV route, with whole-graph propagation,
type-normalized ranks, and the null control, turns the same signature into
specific, robust findings (basal-forebrain anatomy, interneuron-migration
processes). The API route is still useful as a fast first look and as an
independent check against the *current* graph rather than the frozen 2019
snapshot — but interpret its hit counts with the hub bias in mind.

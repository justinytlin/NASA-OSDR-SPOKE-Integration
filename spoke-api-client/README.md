# spoke-api-client

Python client + CLI for the public **SPOKE REST API** (https://spoke.rbvi.ucsf.edu/swagger/).

Unlike the Neo4j/Bolt route (SSH tunnel + credentials), this API is **public, read-only,
and needs no authentication** — it works from anywhere with internet access.

## Requirements

Python 3.10+ with `requests` (`pip install requests`).

## CLI usage

```bash
python3 spoke.py <command> [args] [--format json|csv] [--output FILE]
```

| Command | What it does |
|---|---|
| `version` | SPOKE version + update timestamp for every data source |
| `metagraph` | Node/edge type metagraph (how types interconnect) |
| `types` | Node types, edge types, and available cutoff parameters |
| `search QUERY [--type Disease]` | Lucene name search → `{node_type, identifier, name, score}` |
| `node TYPE ATTR VALUE` | One node's full properties, e.g. `node Disease identifier DOID:2377` |
| `neighborhood TYPE ATTR VALUE` | Neighborhood graph; `--nodes Gene Protein` / `--edges TYPE` filter, `--depth N`, `--cutoff NAME VALUE` |
| `expand TYPE NODE_ID` | Expand a node by internal SPOKE id (from a prior graph result) |
| `sea SMILES_OR_ZINC` | SEA compound-similarity search |

Examples:

```bash
# Find a disease node
python3 spoke.py search "multiple sclerosis" --type Disease

# All genes around it, as CSV
python3 spoke.py --format csv --output ms_genes \
  neighborhood Disease identifier DOID:2377 --nodes Gene

# Tighten the disease-gene text-mining cutoff
python3 spoke.py neighborhood Disease identifier DOID:2377 \
  --nodes Gene --cutoff cutoff_DaG_textmining 5
```

## Library usage

```python
from spoke_api import SpokeClient, split_graph

client = SpokeClient()
hits = client.search("multiple sclerosis", node_type="Disease")
graph = client.neighborhood("Disease", "identifier", hits[0]["identifier"],
                            node_filters=["Gene"])
nodes, edges = split_graph(graph)
```

Graph responses are lists of cytoscape.js-style elements. Node payloads look like
`{"id": ..., "neo4j_type": "Gene", "properties": {...}}`; edge payloads add
`source`/`target` (internal node ids) and `neo4j_type` (e.g. `ASSOCIATES_DaG`).

## Example workflow

[examples/disease_gene_table.py](examples/disease_gene_table.py) searches a disease,
pulls its gene neighborhood, and writes a `disease,gene` association table — the shape
of output you'd join against an external dataset (e.g. differentially expressed genes
from a NASA OSDR study):

```bash
python3 examples/disease_gene_table.py "multiple sclerosis"
```

## Notes

- The API is read-only; integrating external datasets means pulling SPOKE
  subgraphs out and joining them locally.
- Key gene identifiers exposed: Entrez id (`identifier`), symbol (`name`),
  and `ensembl` id — all useful join keys for expression datasets.
- Cutoff parameters (see `types`) control edge stringency, e.g.
  `cutoff_PiP_confidence` (default 0.7), `cutoff_DaG_textmining` (default 3.0).

## examples/ — analyses on the live graph

The API only serves neighborhoods, so these build small subgraphs around a
gene signature instead of using the whole graph. Useful as an independent
check on the PSEV route (`../replication/`) and for anything that needs
current graph content.

| Script | What it does |
|---|---|
| `examples/osdr_spoke_projection.py` | Simplest analysis: signature genes → 1-hop neighborhoods → which diseases/pathways get hit. Fast but hub-biased. |
| `examples/fetch_subgraph.py` | Unions the 1-hop neighborhoods of the signature genes into a local subgraph (`subgraph.jsonl`). |
| `examples/spoke_embeddings_ml.py` | PSEV-style per-sample embeddings on that subgraph (personalized PageRank on up/down-regulated genes) + leave-one-out classification. |
| `examples/welch_nodes.py` | Which nodes *separate* two sample groups: Welch's t-test per node on the embeddings, rank-by-type, top-2.5% flags. |
| `examples/predict_treatments.py`, `alzheimers_multidb.py`, `disease_gene_table.py` | Drug-repurposing-style queries on the live graph. |

```bash
python3 examples/fetch_subgraph.py ...        # genes -> subgraph.jsonl
python3 examples/spoke_embeddings_ml.py ...   # per-sample embeddings
python3 examples/welch_nodes.py --dir OUTDIR
```

Treat `welch_nodes.py` p-values as rankings, not calibrated significance
(embedding dimensions are correlated and sample counts are small). Include
`Anatomy` in the node filters when fetching subgraphs if you want the
anatomy sanity check to work on this route.

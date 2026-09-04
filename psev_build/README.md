# psev_build — gene PSEVs from a current SPOKE Neo4j instance

Rebuilds the Zenodo `gene_psev.zip` (SPOKE v2, 2019; 389,297 nodes, 19,567 gene
PSEVs) against whatever SPOKE version is reachable over Bolt, using the exact
algorithm of the original generator
([BaranziniLab/PSEV](https://github.com/BaranziniLab/PSEV)
`make_psevs_by_node_type.py`, cloned at `../PSEV/`):

| step | original | here |
|---|---|---|
| adjacency | binary, symmetric, from `neo4j_edges.tsv` (both directions present) | binary, symmetric, parallel edges collapsed, self loops dropped |
| transition | column-normalised then transposed = degree-normalised random walk | same, as a sparse CSR |
| restart | one-hot on the gene, jump probability `b = 0.1` | same |
| iteration | power iteration from uniform; stop at L1 change ≤ `a = 0.001` or 40 iterations; renormalise | same defaults; `--converge` runs to 1e-12 |
| output | `gene_group_<g>.tsv` + `raw_psev_0_1_gene_group_<g>_sparse.npy` (1000 genes/group, float64) | same layout; computed in float64, stored as `--dtype float32` by default (half the disk) |

`replication/psev_pipeline.py` (and everything built on `compute_rank_table`)
now accepts a **directory** in this layout wherever it took the zip, so a new
build plugs into the existing pipeline with `--zip <build dir> --spoke <build dir>`.

Validated on the SPOKE v1 graph that ships with the PSEV repo: rows match a
literal port of the original math to ~1e-9 (float32 transition matrix), sum to 1, and put ≈0.10 on the
restart gene (as the Zenodo v2 vectors do).

## Credentials

`export_spoke.py` never prompts. Put these in `psev_build/.env` (or the project
root, or export them):

```
NEO4J_URI=bolt://localhost:7687
NEO4J_USER=neo4j
NEO4J_PASSWORD=...
NEO4J_DATABASE=spoke        # auto-detected from SHOW DATABASES if the name is wrong
```

The instance behind the SSH tunnel (checked 2026-09-02) is Neo4j 2026.07.1,
database `spoke`, a SPOKE snapshot whose newest source timestamp is
2025-03-02: 42.9M nodes / 176M edges, of which Protein = 39.4M (all
organisms; 202k human, 20.4k human Swiss-Prot) and Compound = 2.75M
(381k touch anything other than Compound/Reaction). Four relationship types
hold 93% of the edges and are excluded from the export by default:
`ENCODES_OeP` (Organism→Protein, 38.7M), `PARTOF_PDpP` (51.4M),
`INTERACTS_PiC` (69.2M), `HAS_PhEC` (5.4M). `MiRNA` is exported but worth
excluding from the walk: 2,656 nodes with 1.6M `TARGETS_MtG` edges.

## Run order

```bash
cd psev_build

# 1. how big is it?  (count store, fast)
python3 export_spoke.py stats                       # -> spoke_current/graph_stats.json

# 2. node table (Node = identifier, Node_Name, Node_Type, Node_Index)
python3 export_spoke.py nodes --protein human-reviewed \
        --exclude-types Location,Organism,Version,DatabaseTimestamp,SARSCov2,CellLine,AnatomyCellType,Nutrient

# 3. edge lists, one .npy per relationship type (resumable)
python3 export_spoke.py edges --exclude-rel-types ENCODES_OeP,PARTOF_PDpP,INTERACTS_PiC,HAS_PhEC

# 4. size check, then build
python3 make_gene_psevs.py --spoke spoke_current --dry-run
python3 make_gene_psevs.py --spoke spoke_current --out spoke_current/gene_psev \
        --node-types Gene,Protein,Compound,Disease,Symptom,Anatomy,SideEffect,Pathway,PwGroup,BiologicalProcess,MolecularFunction,CellularComponent,PharmacologicClass,CellType,Complex,Reaction,EC,ProteinDomain,ProteinFamily \
        --prune-types Compound --drop-isolated [--workers 7]
```

`--prune-types Compound` keeps a Compound only if it has an edge to a node of
another type in the universe (drops ChEBI is-a-only and orphan compounds);
`--drop-isolated` removes nodes of any type left without edges. The build
writes its own `node_info.tsv` for the pruned universe, so always pass the
build directory as `--spoke` downstream.

```bash

# 5. integrity (shapes, row sums, restart mass), then same genes new vs 2019
python3 check_build.py spoke_current/gene_psev
python3 compare_psevs.py --a spoke_current/gene_psev --a-nodes spoke_current/gene_psev/node_info.tsv \
        --b ../psev_data/gene_psev.zip --b-nodes ../SPOKE_NASA/spoke_v_2/node_info.tsv --n 30

# 6. use it
python3 ../replication/osdr_to_psev.py TABLE.csv --name OSD-xxx \
        --zip spoke_current/gene_psev --spoke spoke_current/gene_psev
```

## Sizing

Output is `n_genes × n_nodes × 4 bytes` (float32). The 2019 build was
19,567 × 389,297 × 8 B = 61 GB. Published SPOKE sizes since then are
27M nodes / 53M edges (2023 paper) and 42M / 160M (later), which would be
terabytes at full width — `--node-types` restricts the node universe (edges to
dropped types are removed before the walk), which is the intended knob. Compute
scales with `edges × genes`; the v1 test ran 48 genes over 2.1M edges in ~3 s on
4 workers.

## Build of 2026-09-03

`/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03` (external exFAT drive):
19,507 gene PSEVs × 547,873 nodes, float32, 40 GB, 27 min on 7 workers.
Universe = the 19 types in the command above with `--prune-types Compound
--drop-isolated` (2.37M orphan/ChEBI-only compounds and 17.5k isolated nodes
dropped). Downstream outputs live beside it: `genelab_runs_spoke2025-03/`
(six GeneLab rank tables), `meta_out_spoke2025-03_{default,matched}/`,
`osdr_runs_spoke2025-03/` (OSD-564). Logs: `build.log`, `rebuild_and_run.log`.

Gotcha: macOS writes `._*` AppleDouble sidecars next to every file on exFAT;
`PsevStore` ignores them, but delete them (`find DIR -name '._*' -delete`)
before copying the directory anywhere else.

### Results of the 2026-09-03 build vs the 2019 zip

Same gene inputs (six GeneLab studies: 16,172 / 6,698 / 1,599 / 2,642 /
8,275 / 3,698 genes vs 16,181 / 6,701 / 1,598 / 2,643 / 8,276 / 3,701 in 2019),
same pipeline, same Welch meta-analysis. Only 81,553 nodes share an identifier
across the two universes (compounds moved to InChIKeys, proteins were
re-scoped), so every comparison below is over shared nodes.

| | 2019 build | 2025 build |
|---|---|---|
| Gene PSEV vectors, Spearman over shared nodes (20 genes) | | ρ = 0.77–0.88, median 0.87 |
| Welch −log10 p, Space v Ground, Spearman across shared nodes | | ρ = 0.11 (Space v Basal: 0.26) |
| Top-2.5% overlap vs chance: Gene / BiologicalProcess / Protein | | 177 vs 31 / 83 vs 17 / 130 vs 33 |
| Top-2.5% overlap vs chance: Disease / Anatomy / SideEffect / Symptom | | 15 vs 14 / 24 vs 17 / 0 vs 6 / 0 vs 0.5 |
| Paper's Fig-4 showcase nodes at p < 0.025 (7 findable by name) | 7/7 | 5/7 after GO/WikiPathways renames (lost: sympathetic nervous system, taste receptor complex) |
| OSD-564 (right hippocampus) top Anatomy, Space v Ground | diagonal band of Broca, pes | olfactory bulb layers, fusiform gyrus, temporal cortex |

Reading: gene-level and process-level signals carry over between graph
versions well above chance; disease, side-effect and symptom rankings do not,
because those layers were rebuilt almost entirely (Disease 9.1k→11.8k nodes,
Symptom 369→1,405 with new HPO gene edges, SIDER remapped). Treat 2019-build
disease/symptom findings and 2025-build findings as separate analyses rather
than expecting one to confirm the other, and keep the ground-control null
(`replication/osdr_null_control.py`) in the loop for either.

Data quirks in this SPOKE snapshot worth knowing: WikiPathways pathways
appear 2–3× (`WP1531`, `WP1531_r118354`, … are revision copies with identical
gene sets), and raw PSEVs are dominated by Bgee anatomy hubs (blood,
gastrocnemius, heart left ventricle) instead of 2019's brain/testis; the
pipeline's cross-gene z-scoring handles the second, nothing handles the first.

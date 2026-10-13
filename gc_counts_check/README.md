# GC Counts Check — does a control transcriptome point at its own tissue in SPOKE?

Each OSDR dataset is judged **on its own**: its control-sample gene counts are projected onto the gene PSEVs of
the SPOKE 2025-03 snapshot built with only the up-regulation Anatomy–Gene edges, the projection is corrected
for the graph's own baseline, and the resulting Anatomy ranking is checked for the known source tissue. Nothing
is compared across datasets. Run on every mouse bulk RNA-seq control set on OSDR (94 studies, 104 control sets,
31 tissues, 916 samples), the tissue's own gene-linked Anatomy node ranks **first in 50 datasets and in the top
seven in 32 more**; 20 datasets belong to tissues without a gene-linked Anatomy node and are covered by the
cell-type route; two fail. The same pipeline, packaged for any counts table, is
[SPOKE-GeneCounts-Anatomy-Pipeline](https://github.com/justinytlin/SPOKE-GeneCounts-Anatomy-Pipeline).
`METHODS.md` explains every step; `PRESENTER_GUIDE.md` walks through `results/SPOKE_tissue_validation.pptx`.

## Folder layout

```
scripts/gc_counts_to_psev.py    counts -> per-sample SPOKE scores (--manifest, --zip/--spoke build dir, --compact)
scripts/check_linked.py         primary test: rank of the tissue's own node among the 60 gene-linked Anatomy nodes
scripts/check_anatomy.py        secondary test: UBERON family vs all Anatomy nodes (QQ, MWU, KS, per sample)
scripts/celltype_check.py       fallback: rank of the tissue's cell type among the 2,745 CellType nodes
scripts/ontology_families.py    UBERON families (root term + is_a/part_of descendants); families.py = old name regexes
scripts/check_build.py          run the checks for every manifest entry; compile_build.py -> results tables + overview
scripts/plot_celltype.py        cell-type route figures; pilot_stage2.py / plot_stage2.py = edge-set comparison
scripts/survey_osdr.py, build_manifest.py   OSDR survey -> data/_catalog/manifest*.json (which samples are controls)
scripts/inspect_edges.py, inspect_aug_coverage.py   read-only Neo4j queries behind the edge table below
scripts/extract_anatomy_psev.py, gene_hub_rank.py, pilot_stage1.py, check_manifest.py, compile_results.py
                                earlier tooling (PSEV-based linked-node detection, scoring-variant pilot, 2019 run)
run_upauG_all.sh                the full 104-dataset run on the upregulates-only build (three batches)
deck/build.js                   source of results/SPOKE_tissue_validation.pptx (pptxgenjs)
data/                           GeneLab counts + ISA metadata (fetch_data.sh), _catalog/ (manifests), _ontology/, _genesets/
psev_out/                       <name>__2025upAuG_GC_ranks.npz (per-sample scores over all 547,869 nodes, uniform profile,
                                mean weights) and _GC_weights.tsv; not in git
datasets/<name>__2025upAuG/     corrected_linked_2025upAuG60/ (primary), corrected_ontology/ (secondary),
                                corrected_linked_2025upAuG[72]/ (44- and 72-node sets, for reference)
results/                        compiled tables and overview figures (see below); archive_2019/ = results on the 2019 PSEVs
archive/                        earlier cross-dataset analysis; not part of the check
```

## Running the check

The build (40 GB) lives on the external drive: `/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only`
(how it was made: `psev_build/README.md`; only `UPREGULATES_AuG` connects Anatomy to Gene nodes, everything else
in SPOKE is unchanged).

```bash
B=/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only
# one or more manifest entries, one pass over the build (~8 min for one dataset, ~55 min for 13)
python3 scripts/gc_counts_to_psev.py --manifest data/_catalog/manifest_2025upAuG.json --names OSD-245__2025upAuG --compact --zip $B --spoke $B --threads 4
# primary test (60-node set), secondary test, cell-type route
python3 scripts/check_linked.py  --npz psev_out/OSD-245__2025upAuG_GC_ranks.npz --tissue liver --name OSD-245__2025upAuG --spoke $B --linked-ids results/aug_linked_nodes_ge200.tsv --tag 2025upAuG60
python3 scripts/check_anatomy.py --npz psev_out/OSD-245__2025upAuG_GC_ranks.npz --tissue liver --name OSD-245__2025upAuG --spoke $B --corrected --families ontology
python3 scripts/celltype_check.py --tag 2025upAuG --build $B      # all manifest entries -> results/celltype_check_2025upAuG.tsv
# everything at once for the whole manifest, then the tables
python3 scripts/check_build.py --tag 2025upAuG --build $B [--names a,b] [--force]
python3 scripts/compile_build.py --tag 2025upAuG --linked-tag 2025upAuG60
```

For a counts file that is not on OSDR, or for a one-command run, use the standalone repository:
`python run_check.py --counts COUNTS.csv --tissue liver --samples _GC_ --build $B`.

## Datasets and samples

`scripts/survey_osdr.py` lists every *Mus musculus* RNA-seq study on OSDR with a GeneLab
`RSEM_Unnormalized_Counts` file; `scripts/build_manifest.py` picks the control samples from the ISA tables
(spaceflight studies: `Spaceflight` factor exactly "Ground Control"; ground-analog studies: the baseline level
of every condition factor, labelled `analog control`; wild type only; one entry per study × tissue; skipped
studies and reasons in `data/_catalog/manifest_skipped.tsv`). `data/_catalog/manifest_2025upAuG.json` holds
the 104 control sets with explicit sample columns: 83 spaceflight ground-control sets and 21 analog baselines.

| tissue | sets | tissue | sets | tissue | sets |
|---|---|---|---|---|---|
| liver | 14 | adrenal gland | 4 | bone marrow | 2 |
| spleen | 8 | quadriceps | 4 | colon | 2 |
| skin | 7 | lung | 3 | hippocampus | 2 |
| soleus | 7 | cerebellum | 3 | tibialis anterior | 2 |
| kidney | 6 | cerebrum | 3 | extensor digitorum longus | 2 |
| heart | 6 | brain | 2 | eye | 2 |
| gastrocnemius | 5 | | | | |
| retina | 5 | one set each: blood, white adipose tissue, brown adipose tissue, mammary gland, placenta, cecum, bone, mandible, temporal bone, optic nerve | | | |
| thymus | 5 | | | | |

## Method in brief

1. **Weights.** Per sample, CPM from the sample's own library size; genes with CPM ≥ 1 in at least one
   control sample; weight log2(CPM + 1). Mouse Ensembl → mouse Entrez (GeneLab annotation) → human Entrez
   (HomoloGene 68) → SPOKE gene node; 10–14k genes per dataset reach a PSEV.
2. **Projection.** For each gene, z-score its PSEV over the dataset's genes and rank the 547,869 nodes
   within the gene; each sample's score for a node is the expression-weighted mean rank over genes.
3. **Baseline correction.** The same projection with every gene weighted 1 (the uniform profile) is
   subtracted, scaled by the sample's mean weight. What remains is the covariance across genes between
   expression and the node's gene links. Without it the raw ranking is nearly the same in every tissue:
   it correlates with the uniform profile at Spearman 0.93, the raw rankings of liver, kidney, spleen,
   retina, heart, cerebellum and muscle datasets agree at a median of 0.87 (one shared component carries
   71 % of their variance), and the same five kidney-cortex nodes head every one of them. In OSD-245 liver
   the correction moves `right lobe of liver` from 161st to 9th of 13,895 Anatomy nodes (4th to 1st among
   the 60 linked nodes) and puts the five liver-lobule zones at the top; `cortex of kidney` drops from 1st
   to 2,882nd.
4. **Three tests, each within the dataset.** Primary: rank of the tissue's own node among the 60 gene-linked
   Anatomy nodes (CLEAR = 1st, PARTIAL = top 7, NONE below; proxy when only a linked ancestor exists).
   Secondary: the tissue's UBERON family against all Anatomy nodes (QQ plot with a random-node envelope,
   one-sided Mann-Whitney U and KS, per-sample consistency, rank among 14 organ families). Fallback for
   tissues without a linked Anatomy node: rank of the tissue's cell type among the 2,745 CellType nodes.

## Why only the up-regulation edges

Read from the live graph (`scripts/inspect_edges.py`, `inspect_aug_coverage.py`,
`results/anatomy_gene_edge_counts.tsv`):

| edge type | edges | Anatomy nodes | genes | properties |
|---|---|---|---|---|
| `EXPRESSES_AeG` | 1,582,757 | 92 | 19,072 | `sources` only (Bgee); blanket: blood is linked to 18,879 of 19,072 genes |
| `UPREGULATES_AuG` | 85,530 | 72 | 17,018 | `log2fc` (min 1.0, median 1.6, 90th pct 3.3), `pvalue`, `fdr` |
| `DOWNREGULATES_AdG` | 493,600 | 73 | 19,014 | `log2fc`, `pvalue`, `fdr` |

The PSEV walk is binary, so an `expresses` or `downregulates` edge counts as much as an `upregulates` edge,
and the blanket edges turn every linked organ node into a generic hub. Ten pilot datasets, one per tissue,
projected on four builds of the same 2025 graph (`scripts/pilot_stage2.py`, `results/pilot_stage2_<tag>.tsv`,
`results/pilot_stage2_builds_2025only.png`); rank of the tissue's own gene-linked node, number of linked
nodes in brackets:

| tissue | all edges (92 linked) | no 'expresses' (73) | housekeeping removed, all edges (92) | 'upregulates' only (44) |
|---|---|---|---|---|
| retina | no node | no node | no node | no node |
| kidney | 21 | 38 | 20 | no node |
| heart | 8 | 45 | 7 | no node |
| spleen | 75 | 48 | 74 | **2** |
| liver | 35 | 45 | 31 | **1** |
| adrenal gland | 5 | 38 | 5 | **2** |
| lung | 11 | 30 | 11 | **3** |
| cerebellum | 1 | 20 | 1 | **1** |
| white adipose tissue | 10 | 37 | 10 | **1** |
| blood | 82 | 73 | 79 | **1** |
| Spearman(link weight, own rank) | −0.95 | −0.35 | −0.95 | −0.18 |

The 'all edges' and 'upregulates only' columns are computed with the corrected weights; the liver, lung and heart cells of the two middle columns come from the run before the weight fix (the other seven pilot tissues were not affected by it). The pilot's 'upregulates only' column uses the 44 nodes the
PSEV-threshold detector finds; the main results below use the edge-derived 60-node set, on which kidney and heart
become testable.


- **Full graph:** link count predicts own-tissue rank (Spearman −0.95); the most lightly linked nodes top
  almost every dataset whatever the tissue, the most heavily linked ones sink.
- **Without `expresses` but with `downregulates`:** fails everywhere; `putamen` and `amygdala` top every dataset.
- **Housekeeping genes removed (full graph):** reproduces the full-graph column within a few places, so the
  problem is not the gene weighting.
- **`upregulates` only:** every testable tissue's own node in the top three, the top node differs by tissue,
  and the connectivity correlation shrinks to −0.18. This is the configuration used for everything below.

**Which Anatomy nodes are tested.** 72 nodes carry up-regulation edges; the sparsely linked ones (heart 38
edges, brain 19, ectocervix 27) top many datasets because a covariance over a handful of genes is noisy, so
the primary test uses the **60 nodes with ≥ 200 up-regulation edges** (`results/aug_linked_nodes_ge200.tsv`).
For reference, the 44-node set found by the PSEV-threshold detector gives 51 CLEAR / 18 PARTIAL / 35 no node
(`compiled_2025upAuG.tsv`) and the full 72-node set 45 / 27 / 13 NONE / 19 (`compiled_2025upAuG72.tsv`).
Scoring variants tried before the edge set was changed (z-score weighting, a hub-matched permutation null,
housekeeping removal) did not move the failing tissues; the rank-based score was kept (METHODS.md § 15).

## Results: 104 datasets, 60-node set

`results/compiled_2025upAuG60.tsv` (one row per dataset), `results/tissue_summary_2025upAuG60.tsv`,
`results/tissue_overview_2025upAuG60.png`.

| own node | datasets |
|---|---|
| ranked 1st (CLEAR) | 50 |
| ranked 2nd–7th (PARTIAL) | 32 |
| below 7th (NONE) | 2 |
| tissue has no linked Anatomy node | 20 |

| tissue (sets) | own node | 1st | 2nd–7th | median rank | family verdicts (median pct) |
|---|---|---|---|---|---|
| liver (14) | right lobe of liver | 14 | 0 | 1 | CLEAR 14 (0.3) |
| soleus (7), gastrocnemius (5), quadriceps (4), tibialis anterior (2), extensor digitorum longus (2) | gastrocnemius | 20 | 0 | 1 | PARTIAL 20 (5) |
| cerebrum (3), brain (2) | dorsolateral prefrontal cortex | 5 | 0 | 1 | CLEAR 1, PARTIAL 4 (5–7) |
| white adipose tissue (1) | subcutaneous adipose tissue | 1 | 0 | 1 | PARTIAL (31) |
| blood (1) | blood | 1 | 0 | 1 | PARTIAL (0.4) |
| kidney (6) | adult mammalian kidney | 4 | 2 (2nd, behind right lobe of liver) | 1 | PARTIAL 6 (37) |
| cerebellum (3) | cerebellar hemisphere | 2 | 1 | 1 | CLEAR 3 (2.6) |
| spleen (8) | spleen | 2 | 6 (2nd, behind blood) | 2 | PARTIAL 8 (26) |
| skin (7) | skin of abdomen / skin of leg | 1 | 6 (2nd–4th, behind descending thoracic aorta) | 2 | CLEAR 3, PARTIAL 4 (27) |
| heart (6) | heart left ventricle | 0 | 6 (2nd, behind gastrocnemius) | 2 | CLEAR 6 (2.2) |
| adrenal gland (4) | right adrenal gland | 0 | 4 (2nd–3rd, behind descending thoracic aorta) | 2 | PARTIAL 3, NONE 1 (33) |
| hippocampus (2) | Ammon's horn | 0 | 2 (4th, behind prefrontal and cingulate cortex) | 4 | CLEAR 2 (4.7) |
| lung (3) | upper lobe of left lung | 0 | 3 (4th–6th, behind aorta and arteries) | 5 | NONE 3 (68) |
| mammary gland (1) | thoracic mammary gland | 0 | 1 (5th) | 5 | NONE (62) |
| colon (2) | muscle layer of sigmoid colon | 0 | 1 (6th); the other 17th | 11.5 | CLEAR 2 (20) |
| brown adipose tissue (1) | subcutaneous adipose tissue | 0 | 0 (12th) | 12 | NONE (47) |
| retina (5), thymus (5), eye (2), bone marrow (2), placenta, cecum, bone, mandible, temporal bone, optic nerve | no linked node | – | – | – | thymus NONE 5 (84); retina PARTIAL 5 (17); cecum, mandible CLEAR |

- Ground-analog baselines behave like spaceflight ground controls (11 CLEAR / 6 PARTIAL / 4 no node vs
  39 / 26 / 2 NONE / 16).
- Near-misses are a biologically related node winning, and they are consistent: spleen behind `blood` in
  6 of 8, heart behind `gastrocnemius` in 6 of 6, adrenal gland and skin behind `descending thoracic aorta`.
- Residual hubs: `gastrocnemius` is in the top four of 43 datasets, `heart left ventricle` of 29,
  `descending thoracic aorta` of 28.
- The secondary (family) test agrees where families are clean (liver, heart, cerebellum, hippocampus) and
  is weaker where the UBERON family mixes parenchyma with vessels and skin appendages (kidney, spleen,
  adrenal gland).

## Cell-type route for tissues without a gene-linked Anatomy node

The graph has 2,745 CellType nodes, 59 with gene edges (33,827 `EXPRESSEDIN_GeiCT`) and 846 `PARTOF_CTpA`
links to Anatomy nodes. `scripts/celltype_check.py --tag 2025upAuG --build $B` ranks them per sample by the
same corrected score (`results/celltype_check_2025upAuG.tsv`); `scripts/plot_celltype.py --tag 2025upAuG`
draws the 20 datasets without a linked Anatomy node (`results/celltype_nolinked_2025upAuG.png`) and, with
`--all`, all 104 (`results/celltype_all_2025upAuG.png`). Ranks are out of 2,745.

| tissue (sets) | best matching cell type | rank, median [range] |
|---|---|---|
| retina (5) | camera-type eye photoreceptor cell | 1 [1–2] |
| bone marrow (2) | CD115-positive monocyte | 1 |
| placenta (1) | syncytiotrophoblast cell | 1 |
| cecum (1) | Paneth cell | 2 |
| thymus (5) | T cell | 17 [12–19] |
| eye (2) | retinal rod cell | 19 [15–23] |
| optic nerve, mandible, bone, temporal bone | (oligodendrocyte / chondrocyte) | 30, 47, 121, 146 |

For tissues that do have an Anatomy node the cell-type ranking agrees with it: hepatocyte 1st in 14 of 14
liver sets, proximal tubular cell 1st in 6 of 6 kidney sets, suprabasal keratinocyte 1st in 4 of 7 skin sets,
colonic smooth muscle fibre 1st in both colon sets, cardiac pacemaker cell 3rd–4th in heart, sinoatrial
myocyte 3rd in the muscles; adrenal gland has no steroidogenic cell type among the 59 and ranks a kidney cell
at ~109. Residual hub: sinoatrial-node myocytes and epicardial adipocytes are in the top five of 27 % of datasets.

## Tissue coverage in three groups

| group | tissues | how |
|---|---|---|
| validated directly (20) | liver, kidney, heart, spleen, lung, adrenal gland, skin, brain, cerebrum, cerebellum, hippocampus, white adipose tissue, blood, mammary gland, the five hindlimb muscles, colon (weak) | own Anatomy node among the 60 linked nodes |
| validated by cell type (4) | retina, bone marrow, placenta, cecum; thymus weak (T cell median 17th), eye marginal | own CellType node among 2,745 |
| not testable (4) | bone, mandible, temporal bone, optic nerve | no Anatomy or CellType node with gene edges |

So 27 of 31 tissues have a direct within-dataset test; thymus, eye, colon and brown adipose tissue carry an
explicit "weak" label. The gaps follow the human atlases behind SPOKE's edges: GTEx never sampled thymus,
retina, bone marrow or placenta, and no atlas edge exists for bone.

## Reproducibility note (2026-10-13)

While porting the pipeline to its standalone repository, the CPM step was found to return wrong values for
frames above one million cells: pandas routes such elementwise operations to numexpr, and numexpr 2.10.1
under numpy 2.4.6 computed them incorrectly with no warning. The 13 control sets with 19 or more samples
(OSD-243, 244, 245, 248, 253, 254, 295, 379, 580, 720, 771, 825, 881) had 13–21 % of their weight cells
wrong; the other 91 were exact (`results/weight_audit_numexpr_2025upAuG.tsv`). The CPM step now uses explicit
numpy arithmetic with the pandas fast path disabled, the 13 datasets were re-projected and re-checked, and
every table and figure in `results/` was regenerated. Two verdicts changed, OSD-243 skin and OSD-771 kidney
(own node 2nd instead of 1st); the other eleven were unchanged.

## Archive

`results/archive_2019/` holds the tables and overview figures of the same datasets on the 2019 Zenodo PSEVs
(the development history that led to the baseline correction and to the edge-set change; superseded, and
partly computed before the fix above). `archive/cross_dataset/` holds the earlier analysis that compared
datasets with one another; `archive/combined_check_figures/` the previous three-panel figures.

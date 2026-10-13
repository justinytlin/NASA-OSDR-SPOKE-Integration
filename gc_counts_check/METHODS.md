# Methods guide: mapping gene counts onto SPOKE and checking tissue identity

What `scripts/gc_counts_to_psev.py`, `scripts/check_linked.py`, `scripts/check_anatomy.py` and
`scripts/celltype_check.py` do, in the order they do it, on the gene PSEVs of the SPOKE 2025-03 snapshot built
with only the up-regulation Anatomy–Gene edges. Every dataset is processed and judged on its own; nothing is
compared across datasets. The same pipeline, packaged for reuse, is
[SPOKE-GeneCounts-Anatomy-Pipeline](https://github.com/justinytlin/SPOKE-GeneCounts-Anatomy-Pipeline).

## 1. Inputs

| input | source | role |
|---|---|---|
| GeneLab RSEM unnormalized counts (`GLDS-xxx_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv`) | OSDR, GeneLab bulk RNA-seq pipeline | rows = mouse Ensembl gene ids, columns = samples, values = expected read counts |
| Gene-PSEV build (`/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only/`) | `psev_build/make_gene_psevs.py` on the SPOKE 2025-03 Neo4j export | `node_info.tsv` (547,869 nodes: 19,507 Gene, 13,895 Anatomy, 2,745 CellType, ...) and one 547,869-long float32 vector per gene in 20 groups |
| GeneLab annotation (ENSEMBL, ENTREZID columns of an OSD-564 differential-expression table) | OSDR | mouse Ensembl → mouse Entrez |
| NCBI HomoloGene build 68 (`psev_data/homologene_build68.data`) | NCBI (frozen copy) | mouse Entrez → human Entrez |
| UBERON (`data/_ontology/uberon-basic.obo`) | OBO Foundry | tissue families (section 9) |
| `results/aug_linked_nodes_ge200.tsv` | live SPOKE graph (section 7) | the 60 Anatomy nodes with ≥ 200 up-regulation edges |

Which samples count as controls for each study is in section 12.

## 2. From counts to per-gene weights

For each selected sample the counts are converted to counts per million using that sample's own library
size (the column sum over all genes in the file), and a gene is kept if it reaches CPM ≥ 1 in at least one
selected sample. The weight of gene *g* in sample *s* is

    w[g, s] = log2(CPM[g, s] + 1)

so highly expressed genes weigh more, unexpressed genes weigh nothing, and the scale is compressed the way
expression data usually is. All of this arithmetic is done on numpy arrays: pandas' elementwise fast path
(numexpr) returned silently wrong values on large frames in the development environment (section 16).

## 3. Mapping mouse genes to SPOKE gene nodes

SPOKE genes are human Entrez ids, so each mouse Ensembl id goes through two lookups:

1. **Ensembl → mouse Entrez** using the GeneLab annotation columns. Multi-mapped entries such as
   `54204|100043580` are split so every Entrez id is tried.
2. **mouse Entrez → human Entrez** through HomoloGene build 68: two genes are orthologs when they share a
   HomoloGene group id.
3. Keep only human ids that exist as Gene nodes of the build.

When several mouse genes map to the same SPOKE gene, their weights are averaged. Typical yield: 55–57k
Ensembl rows → 13–15k expressed genes → 10–14k SPOKE gene nodes with a PSEV.

## 4. What a PSEV is, and which graph it was computed on

A Propagated SPOKE Entry Vector (Nelson et al. 2021) is a random walk with restart over the whole SPOKE
graph, started from one gene node with restart probability 0.1 on a binary, degree-normalised adjacency
matrix. Its entries are the stationary visiting probabilities of all nodes, so each gene's PSEV says how
strongly every node in SPOKE is connected to that gene through the network.

The vectors used here were computed by `psev_build/make_gene_psevs.py` from the 2025-03 SPOKE snapshot with
Compound nodes pruned to those connected to another type, isolated nodes dropped, and two Anatomy–Gene
relationship types left out: `EXPRESSES_AeG` and `DOWNREGULATES_AdG`. Only `UPREGULATES_AuG` connects
Anatomy to Gene nodes in this build. Every other edge type of SPOKE is unchanged. Section 7 explains why.

## 5. Projecting a sample onto SPOKE (the two passes)

For one dataset with gene set *G* (the mapped genes) and samples *S*, the build is read once, group by group:

**Pass 1 (normalisation).** For every node *n*, compute the mean μ[n] and standard deviation σ[n] of
PSEV[g, n] over the genes g in *G*. This removes the "hub" component: nodes that every gene's walk visits
often get a high mean, and only deviations from it count.

**Pass 2 (rank and weight).** For each gene g:

    z[g, n]    = (PSEV[g, n] − μ[n]) / σ[n]
    rank[g, n] = position of node n when the 547,869 z-scores of gene g are sorted ascending
                 (0 = least associated with g, 547,868 = most associated)

then for each sample s

    score[s, n] = Σ_g  w[g, s] · rank[g, n]  /  |G|

so a node scores high in a sample when it is strongly linked to many of that sample's highly expressed
genes. Converting z-scores to within-gene ranks before weighting follows the original pipeline and keeps any
single gene with an extreme PSEV from dominating.

`gc_counts_to_psev.py --compact` carries one extra pseudo-sample whose weights are all 1 (section 6) and
writes `psev_out/<name>_GC_ranks.npz` with `scores` (samples × 547,869, float32), `uniform`, `mean_weight`,
`samples`, `n_genes`, plus `<name>_GC_anatomy_ranks.tsv` with the raw per-type ranks of the Anatomy rows.
All datasets of a manifest are processed in one streaming pass; the per-dataset arithmetic is identical.

## 6. Baseline correction

With all-positive weights, most of score[s, n] does not depend on the tissue. On this build the raw Anatomy
ranking of a liver dataset correlates with the all-genes-weighted-equally profile at Spearman 0.93, and the
raw rankings of liver, kidney, spleen, retina, heart, cerebellum and muscle datasets correlate with one
another at a median of 0.87, with one shared component carrying 71 % of their variance. The same five
kidney-cortex nodes head the raw ranking of every one of those datasets, whatever the tissue. That profile is
set by the shape of each node's PSEV distribution across genes, not by expression.

The tissue-independent part is what the score would be if every gene had the same weight:

    uniform[n] = (1/|G|) Σ_g rank[g, n]

The corrected score subtracts it, scaled by the sample's mean weight:

    corrected[s, n] = score[s, n] − mean_g(w[g, s]) · uniform[n]
                    = (1/|G|) Σ_g (w[g, s] − mean_g w[·, s]) · rank[g, n]

which is the covariance, across the dataset's genes, between how strongly a gene is expressed in the sample
and how strongly that gene is linked to node *n*. It is positive when the sample's highly expressed genes are
the ones SPOKE connects to the node, zero when expression and connectivity are unrelated, and negative when
the node's genes are the weakly expressed ones. It uses only the dataset's own genes and samples.

Example, OSD-245 liver (19 ground controls): raw, `cortex of kidney` ranks 1st of 13,895 Anatomy nodes and
`right lobe of liver` 161st; corrected, the five liver-lobule zones are the top five, `right lobe of liver`
is 9th and `cortex of kidney` 2,882nd. Among the 60 gene-linked nodes `right lobe of liver` moves from 4th
to 1st. Every check below uses the corrected score.

What the correction does not remove: genes that are highly expressed in every tissue also tend to be well
connected in SPOKE, so nodes linked to ubiquitous genes keep a mild advantage in every tissue (the residual
hubs in the README).

## 7. Which Anatomy nodes can be tested, and why only up-regulation edges

SPOKE loads all of UBERON, but Anatomy–Gene edges exist for only a few dozen tissues, taken from human
expression atlases. Read from the live graph (`scripts/inspect_edges.py`, `inspect_aug_coverage.py`,
`results/anatomy_gene_edge_counts.tsv`):

| edge type | edges | Anatomy nodes | genes | properties |
|---|---|---|---|---|
| `EXPRESSES_AeG` | 1,582,757 | 92 | 19,072 | `sources` only; blanket: blood is linked to 18,879 of 19,072 genes |
| `UPREGULATES_AuG` | 85,530 | 72 | 17,018 | `log2fc` (min 1.0, median 1.6), `pvalue`, `fdr` |
| `DOWNREGULATES_AdG` | 493,600 | 73 | 19,014 | `log2fc`, `pvalue`, `fdr` |

The PSEV walk is binary, so an `expresses` or `downregulates` edge counts exactly as much as an
`upregulates` edge. On the full graph the blanket edges turn every linked organ node into a generic hub: a
node's rank in its own tissue is predicted by how many edges it has (Spearman −0.95 between link weight and
own-tissue rank across the pilot tissues), the most lightly linked nodes top almost every dataset, and the
heavily linked ones sink. Removing only the `expresses` edges leaves the `downregulates` edges, which push
the wrong way for an unsigned walk, and fails everywhere. Keeping only the `upregulates` edges removes most of the
artifact (correlation −0.18) and puts every testable tissue's own node in the top three. A rebuild that
filters `expresses` edges by level is not possible: they carry no value. Section 14 points to the numbers.

Every other Anatomy node scores second-hand through `ISA`/`PARTOF` edges to a linked neighbour, so the
primary test ranks the linked nodes only. With all 72 up-regulation-linked nodes, the sparsely linked ones
(heart, 38 edges; brain, 19) top many datasets because a covariance over a handful of genes is noisy; the
set used is the **60 nodes with ≥ 200 up-regulation edges** (`results/aug_linked_nodes_ge200.tsv`).

## 8. Primary test: rank of the tissue's own linked node

`scripts/check_linked.py --npz <name>_GC_ranks.npz --tissue <family> --linked-ids results/aug_linked_nodes_ge200.tsv --tag 2025upAuG60`
ranks the 60 linked nodes by corrected score within each sample (1 = top), averages the ranks over samples,
and reports the rank of the tissue's own linked node(s): the linked nodes inside its UBERON family
(section 9), or, if there are none, the nearest linked ancestors of the family root (flagged proxy).
CLEAR = a target ranked 1st; PARTIAL = best target in the top 7; NONE otherwise; NO LINKED NODE when neither
the family nor any ancestor is linked, in which case the cell-type route (section 11) decides. Figures:
`A_linked_nodes` (top 20 of the 60 with the targets highlighted) and `C_per_sample` (rank of the best
target in every sample). Output: `datasets/<name>__2025upAuG/corrected_linked_2025upAuG60/`.

## 9. Tissue families from UBERON

`scripts/ontology_families.py` reads UBERON and defines each family as a root term (liver, kidney, retina,
hippocampal formation, telencephalon for cerebrum, cranial nerve II for optic nerve, skeletal muscle
tissue + skeletal muscle organ for every hindlimb muscle, ...) plus every term below it through `is_a` and
`part_of`, intersected with the build's Anatomy nodes. `develops_from` is not followed, so primordia are
excluded. The 2025 Anatomy vocabulary with gene links is GTEx-like (`right lobe of liver`, `upper lobe of
left lung`, `skin of leg`) and the families pick those nodes up through `part_of`.

## 10. Secondary test: the whole family against all Anatomy nodes

`scripts/check_anatomy.py --npz <name>_GC_ranks.npz --corrected --families ontology` ranks all 13,895
Anatomy nodes by corrected score within each sample and converts the rank to a percentile
(0 = most associated). The **pooled** percentile of a node is its mean over samples. If nothing about the
sample favoured any tissue the pooled percentiles would be spread uniformly between 0 and 100; that uniform
spread is the null every figure is compared against.

**Figure A, QQ plot.** Take the *k* pooled percentiles of the family's nodes and sort them; plot the *i*-th
smallest against 100·(i − 0.5)/k. On the diagonal the family is indistinguishable from random nodes; below
it the family ranks nearer the top. Thin lines repeat this per sample; gray lines show the other organ
families in the same dataset. The gray band is the pointwise 2.5–97.5 % range of 2,000 random sets of *k*
Anatomy nodes drawn from the dataset's own pooled percentiles.

- *Mann-Whitney U* (`scipy.stats.mannwhitneyu`, alternative `less`): is a randomly chosen family node
  ranked above a randomly chosen other node more often than chance? Driven by the family's median.
- *Kolmogorov-Smirnov* (`scipy.stats.ks_2samp`, alternative `greater`): is the family's cumulative
  distribution above the others'? Sensitive to the largest gap anywhere, so it can be significant when only
  part of the family (the parenchyma, say) is shifted.

Both are one-sided because the question is directional.

**Figure B, organ families.** The family's median percentile against 14 organ families defined the same way
(hippocampus, thymus, spleen, liver, cerebral cortex, heart, kidney, lung, skeletal muscle, bone marrow,
skin, intestine, testis, blood; the dataset's own tissue is added when it is not one of them); rows are
sorted by median, so the source family's position is its "rank of N". This is the within-dataset
replacement for comparing datasets to each other.

**Figure C, per sample.** Each sample's family median and its own Mann-Whitney p; filled points are
significant at 0.05.

**Verdict.** CLEAR = pooled Mann-Whitney p < 0.01, KS p < 0.01, ≥ 80 % of samples individually p < 0.05,
family ranked first among the organ families; PARTIAL = pooled p < 0.05 but a criterion fails; NONE = no
pooled shift. Families mix parenchyma with vessels, nerves and primordia, so this test is reported next to
the primary one rather than driving the dataset verdict. Output: `datasets/<name>__2025upAuG/corrected_ontology/`.

## 11. Cell-type route

The graph has 2,745 CellType nodes, 59 of them with `EXPRESSEDIN_GeiCT` gene edges (33,827 edges, a single
marker-gene source) and 846 `PARTOF_CTpA` links to Anatomy nodes. `scripts/celltype_check.py --tag 2025upAuG
--build <dir>` takes the corrected scores already stored in `<name>_GC_ranks.npz`, keeps the CellType
columns, ranks them per sample (1 = top), pools the ranks, and reports the best-ranked cell type whose name
matches the tissue's pattern (`FAM` in the script) with the same 1 / top-7 thresholds. No projection is
re-run. `scripts/plot_celltype.py` draws the per-dataset figure for the tissues without a linked Anatomy
node. This is what makes retina, bone marrow, placenta and cecum testable.

## 12. Which datasets and which samples

`scripts/survey_osdr.py` queries OSDR for every *Mus musculus* study with an RNA-seq assay, caches each
study's file listing and ISA metadata, and keeps those with a GeneLab `RSEM_Unnormalized_Counts` file.
`scripts/build_manifest.py` turns the sample tables into a manifest with one entry per (study, tissue):

- **Control samples.** For spaceflight studies, samples whose `Spaceflight` factor is exactly "Ground
  Control" (basal, vivarium, cohort and rerun controls are not used). For ground-analog studies (hindlimb
  unloading, radiation, diet), samples at the baseline level of every condition factor; these are labelled
  `analog control`. Genotype is restricted to wild type / control. Age, sex, strain, duration, dissection
  and mission are pooled.
- **Tissue.** From the `Tissue` factor when a study has several tissues, otherwise from
  `Characteristics[Material Type]`; left/right and regional names are merged. Multi-tissue studies are split
  into one entry per tissue (`OSD-457_spleen`).
- **Skipped**, with the reason in `data/_catalog/manifest_skipped.tsv`: cell cultures, tumours, plasma,
  methods studies without a control group, ISA sample names that do not match the counts columns, fewer
  than two control samples.

`data/_catalog/manifest_2025upAuG.json` lists the 104 control sets (90 manifest entries plus the first 14
studies with their original `_GC_` column selection) with explicit sample columns.

## 13. Compiled tables and figures

`scripts/check_build.py --tag 2025upAuG --build <dir>` runs the ontology check for every entry;
`check_linked.py --linked-ids` is run per entry for the 60-node set; `scripts/compile_build.py --tag 2025upAuG
--linked-tag 2025upAuG60` writes `results/compiled_2025upAuG60.tsv` (one row per dataset),
`results/tissue_summary_2025upAuG60.tsv` (one row per tissue) and `results/tissue_overview_2025upAuG60.png`;
`celltype_check.py` writes `results/celltype_check_2025upAuG.tsv` and `plot_celltype.py` the two cell-type
figures. These tabulate per-dataset verdicts; they do not test datasets against each other.

## 14. The edge-set comparison

The same scripts run on any gene-PSEV build directory (`--zip <dir> --spoke <dir>`; the build's own
`node_info.tsv` defines the node universe). Ten pilot datasets, one per tissue (retina, kidney, heart,
spleen, liver, adrenal gland, lung, cerebellum, white adipose tissue, blood), were projected on four 2025
builds: the full graph, the graph without `expresses` edges, the graph with only `upregulates` edges, and
the full graph with housekeeping genes removed from the counts (HRT Atlas mouse list,
`data/_genesets/hk_mouse_ensembl.txt`). `pilot_stage2.py --tag <tag> --build <dir> --run-checks` writes
`results/pilot_stage2_<tag>.tsv`; `plot_stage2.py --no-2019` draws `results/pilot_stage2_builds_2025only.png`.
The README carries the table.

## 15. Scoring variants that were tested and rejected

Before the edge set was changed, three ways of altering the scoring were tried on the same pilot tissues,
each dataset alone and baseline-corrected: weighting the genes' z-scores directly instead of their
within-gene ranks (`--zscore-arm`, `--score-key zscore`); converting the corrected score to a z against
50 permutations of each sample's weights among genes of similar hub rank (`--perm-cov`, `--score-key perm`);
and removing housekeeping genes before weighting (`exclude_genes` in the manifest). None moved the failing
tissues' own nodes by more than a few places, and housekeeping removal reproduced the unmodified method on
the 2025 graph as well (`2025nohk` in section 14). The rank-based score with baseline correction was kept.

## 16. Reproducibility note: the pandas/numexpr fault and the re-run

While porting the pipeline to its standalone repository (2026-10-13) the CPM step was found to return
wrong values for frames above one million cells: pandas routes such elementwise operations to numexpr, and
numexpr 2.10.1 under numpy 2.4.6 computed them incorrectly without any warning. In this project that
affected the 13 control sets with 19 or more samples (OSD-243, 244, 245, 248, 253, 254, 295, 379, 580,
720, 771, 825, 881), where 13–21 % of the weight cells were wrong; the other 91 datasets were exact
(`results/weight_audit_numexpr_2025upAuG.tsv`). The CPM step now uses explicit numpy arithmetic and the
pandas fast path is disabled; the 13 datasets were re-projected and re-checked on the upregulates-only
build and every compiled table and figure in `results/` was regenerated. Two verdicts changed
(OSD-243 skin and OSD-771 kidney, CLEAR → PARTIAL); the rest were unchanged.

## 17. Reading the results

- **Near-misses are usually a related node winning**: spleen behind blood, heart behind gastrocnemius,
  lung behind an artery. `report`/`summary.json` list the top of each ranking so this is visible.
- **Coverage follows the human atlases.** Thymus, retina, bone marrow and placenta have no Anatomy node
  with gene edges; the cell-type route covers three of them, thymus only weakly. Bone, mandible, temporal
  bone and optic nerve have no node of either kind with gene edges.
- **Residual hubs.** `gastrocnemius` is in the top four of about half the datasets; sinoatrial-node
  myocytes recur in cell-type rankings.
- **Small n.** With two or three control samples the per-sample criteria are fragile.

## References

- Nelson CA et al. (2021) *Life* 11:42 — the PSEV method and the GeneLab-to-SPOKE pipeline whose code this
  follows (`replication/psev_pipeline.py` is the faithful port; `psev_build/` rebuilds the vectors).
- Morris JH et al. (2023) *Bioinformatics* 39:btad080 — SPOKE.
- Mungall CJ et al. (2012) *Genome Biol* 13:R5 — UBERON.
- NCBI HomoloGene build 68 — mouse/human orthology.
- GeneLab bulk RNA-seq processing pipeline (GLbulkRNAseq) — source of the RSEM counts and gene annotation.
- Mann HB, Whitney DR (1947); Kolmogorov (1933) / Smirnov (1939) — the two tests, as implemented in SciPy.

# Methods guide: mapping ground-control gene counts onto SPOKE and checking tissue separation

This describes exactly what `scripts/gc_counts_to_psev.py` and `scripts/check_anatomy.py` do, in
the order they do it. Every dataset is processed and judged on its own; nothing is compared
across datasets.

## 1. Inputs

| input | source | role |
|---|---|---|
| GeneLab RSEM unnormalized counts (`GLDS-xxx_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv`) | OSDR, GeneLab bulk RNA-seq pipeline | rows = mouse Ensembl gene ids, columns = samples, values = expected read counts |
| SPOKE v2 node table (`SPOKE_NASA/spoke_v_2/node_info.tsv`) | Nelson et al. 2021 code repository | the 389,297 nodes of the 2019 SPOKE (19,567 Gene nodes, 13,257 Anatomy nodes, ...) |
| Gene PSEVs (`psev_data/gene_psev.zip`, 25.6 GB) | Zenodo 10.5281/zenodo.4404618 | one 389,297-long vector per SPOKE gene node |
| GeneLab annotation (ENSEMBL, ENTREZID columns of an OSD-564 differential-expression table) | OSDR | mouse Ensembl → mouse Entrez |
| NCBI HomoloGene build 68 (`psev_data/homologene_build68.data`) | NCBI (frozen copy) | mouse Entrez → human Entrez |

Only ground-control samples are used: columns whose name contains `_GC_` (for OSD-168, which
mixes RR-1 and RR-3 animals, `_RR3_GC_`).

## 2. From counts to per-gene weights

For each GC sample the counts are converted to counts per million using that sample's own
library size (the column sum over all genes in the file), and a gene is kept if it reaches
CPM ≥ 1 in at least one GC sample. The weight of gene *g* in sample *s* is

    w[g, s] = log2(CPM[g, s] + 1)

so highly expressed genes weigh more, unexpressed genes weigh nothing, and the scale is
compressed the way expression data usually is. This is the "expression-weighted" extension
of Nelson et al. (2021), whose published method weights genes by log2 fold change of a
contrast instead; the arithmetic downstream is identical.

## 3. Mapping mouse genes to SPOKE gene nodes

SPOKE genes are human Entrez ids, so each mouse Ensembl id goes through two lookups:

1. **Ensembl → mouse Entrez** using the GeneLab annotation columns. Multi-mapped entries such
   as `54204|100043580` are split so every Entrez id is tried.
2. **mouse Entrez → human Entrez** through HomoloGene build 68: two genes are orthologs when
   they share a HomoloGene group id. The 2020-era MGI report used by the original paper has
   since dropped its HomoloGene column, so the frozen NCBI build is used instead.
3. Keep only human ids that exist as Gene nodes in SPOKE v2.

When several mouse genes map to the same SPOKE gene, their weights are averaged. Typical
yield: 55–57k Ensembl rows → 13–20k expressed genes → 11–14k SPOKE gene nodes.

## 4. What a PSEV is

A Propagated SPOKE Entry Vector (Nelson et al. 2021) is a random walk with restart over the
whole SPOKE graph, started from one gene node with restart probability 0.1 on a binary,
degree-normalised adjacency matrix. Its entries are the stationary visiting probabilities of
all 389,297 nodes, so each gene's PSEV says how strongly every node in SPOKE is connected to
that gene through the network. The vectors are precomputed and read from the Zenodo archive
in 20 groups.

## 5. Projecting a sample onto SPOKE (the two passes)

For one dataset with gene set *G* (the mapped genes) and samples *S*:

**Pass 1 (normalisation).** For every node *n*, compute the mean μ[n] and standard deviation
σ[n] of PSEV[g, n] over the genes g in *G*. This removes the "hub" component: nodes that every
gene's walk visits often (very connected nodes) get a high mean, and only deviations from it
count.

**Pass 2 (rank and weight).** For each gene g:

    z[g, n]    = (PSEV[g, n] − μ[n]) / σ[n]
    rank[g, n] = position of node n when the 389,297 z-scores of gene g are sorted ascending
                 (0 = least associated with g, 389,296 = most associated)

then for each sample s

    score[s, n] = Σ_g  w[g, s] · rank[g, n]  /  |G|

so a node scores high in a sample when it is strongly linked to many of that sample's highly
expressed genes. Converting z-scores to within-gene ranks before weighting follows the
original pipeline and keeps any single gene with an extreme PSEV from dominating.

**Output ranks.** For every sample the nodes are ranked by score, overall and within each node
type. The per-type rank `Rank_by_type_<sample>` runs from 1 (highest score = most associated
with the sample) to the size of the type (13,257 for Anatomy). The rank table has one row per
SPOKE node and one such column per sample.

All datasets are processed in one streaming pass over the archive; the per-dataset arithmetic
is unchanged, and the result matches the earlier single-dataset run to Spearman 0.999999.

## 6. From ranks to anatomy percentiles

`check_anatomy.py` takes the Anatomy rows only and turns each sample's rank-by-type into a
percentile:

    pct[n, s] = 100 · Rank_by_type[n, s] / 13,257

so 0 means the node is the Anatomy node most associated with the sample's expressed genes and
50 is the middle of the pack. The **pooled** percentile of a node is the mean of pct[n, s]
over the dataset's samples. If nothing about the sample favoured any tissue, these pooled
percentiles would be spread uniformly between 0 and 100 across nodes; that uniform spread is
the null every figure is compared against.

## 7. Tissue families

SPOKE's Anatomy nodes are UBERON terms and a tissue is represented by many of them (the spleen
by "spleen", "white pulp of spleen", "splenic artery", "spleen primordium", ...). A **family**
is the set of Anatomy nodes whose name matches an include regex and not an exclude regex,
defined in `scripts/families.py`. Examples: hippocampus = `hippocamp|dentate gyrus|Ammon|
subiculum` (51 nodes); spleen = `spleen|splenic` minus `flexure` (31); liver = `liver|hepat`
minus hepatopancreas/flexure terms (85); kidney = `kidney|renal|nephron|glomerul` minus adrenal
and olfactory/cerebellar glomerulus terms (167). Families are name-based and therefore include
vessels, ducts and embryonic primordia alongside parenchyma; this matters for interpretation
(section 12).

## 8. Figure A: the QQ plot

**Question.** Do the family's nodes sit higher in the ranking than an equally sized random set
of Anatomy nodes would?

**Construction.** Take the *k* pooled percentiles of the family's nodes and sort them. Plot the
*i*-th smallest against its expected value under the uniform null, 100·(i − 0.5)/k. If the
family were indistinguishable from random nodes the points would lie on the diagonal. A curve
below the diagonal means the family's nodes have smaller percentiles than expected, i.e. they
rank nearer the top. Thin lines repeat this for each individual sample; gray lines show the
other organ families in the same dataset for scale.

**Null envelope.** 2,000 sets of *k* nodes are drawn at random (without replacement) from all
13,257 Anatomy nodes, each set is sorted, and the pointwise 2.5th and 97.5th percentiles of
those sorted draws form the gray band. Because it is built from the dataset's own pooled
percentiles it reflects whatever irregularities that ranking has. A curve that leaves the band
is beyond what chance produces for a family of that size.

**Tests.**
- *Mann-Whitney U* (`scipy.stats.mannwhitneyu`, alternative `less`): compares the family's
  percentiles with those of all other Anatomy nodes; the alternative is that a randomly chosen
  family node has a smaller percentile than a randomly chosen non-family node. It is a test on
  the whole distribution's location, so it is driven by the family's median.
- *Kolmogorov-Smirnov* (`scipy.stats.ks_2samp`, alternative `greater`): compares the two
  empirical cumulative distributions; the alternative is that the family's CDF lies above the
  others', which is the case when family percentiles are stochastically smaller. KS is
  sensitive to the largest gap between the curves wherever it occurs, so it can be significant
  when only part of the family (say the parenchymal terms) is shifted while the median is not.

Both are one-sided because the question is directional: does the tissue rank *higher*.

## 9. Figure B: organ families

**Question.** Within this dataset, is the source tissue the family that ranks best?

**Construction.** For each of 13 fixed organ families (hippocampus, thymus, spleen, liver,
cerebral cortex, heart, kidney, lung, skeletal muscle, bone marrow, skin, intestine, testis,
blood, with the dataset's own tissue added if it is not one of them) every node is drawn as a
dot at its pooled percentile, and a black bar marks the family's median. Families are sorted by
median, so the row order is the ranking and the source family's position is its "rank of N".
The dashed line at 50 marks no association. The same Mann-Whitney test as in Figure A is run for
every family and stored in `organ_families.tsv`.

This panel is the within-dataset replacement for comparing datasets to each other: instead of
asking whether liver scores better in a liver dataset than in a spleen dataset, it asks
whether liver scores better than spleen inside the liver dataset.

## 10. Figure C: per-sample consistency

**Question.** Is the pooled result carried by every sample or by a few?

For each sample separately, the median percentile of the family's nodes is plotted, and the
Mann-Whitney test (family vs other Anatomy nodes, alternative `less`) is run on that sample's
own percentiles. Filled points are samples significant at 0.05; the horizontal line is the
pooled median. `per_sample.tsv` holds the numbers.

## 11. The verdict

`summary.json` records a verdict from four criteria (thresholds are command-line flags):

| criterion | default |
|---|---|
| pooled Mann-Whitney p (family vs other Anatomy nodes) | < 0.01 |
| pooled KS p | < 0.01 |
| fraction of samples individually significant at 0.05 | ≥ 0.8 |
| family rank among the organ families | 1 |

CLEAR = all four pass. PARTIAL = pooled Mann-Whitney p < 0.05 but at least one criterion fails.
NONE = no pooled shift. The JSON also lists the family's best-ranked nodes and, for context, the
top-ranked real Anatomy nodes after removing UBERON life-stage and organism pseudo-nodes.

## 12. Reading the results

- **Absolute position is not the whole story.** On the 2019 build the very top of every Anatomy
  ranking is occupied by life-stage pseudo-nodes and brain/neural hub nodes whatever the tissue,
  so the check asks where the tissue family sits relative to chance and to other organs, not
  whether it is at rank 1.
- **Family definition drives PARTIAL/NONE.** In nearly every dataset the tissue's parenchymal
  terms are in the top 1–3 % (liver lobule, zona fasciculata, renal papilla, corneal epithelium,
  Purkinje fiber) while vessels, ducts and primordia in the same name-based family sit at the
  bottom. A narrower family (for example `--tissue "liver parenchyma"`, or `--family-regex`)
  tests the parenchyma alone; the broad family tests the whole UBERON neighbourhood.
- **Small n.** With two or three GC samples the per-sample criterion is fragile.
- **Node-level anti-association is informative too.** Retina nodes rank at the bottom in retina
  samples; a family that lands consistently below the diagonal is telling you about the graph's
  wiring, not the sample.

## References

- Nelson CA et al. (2021) *Life* 11:42 — the PSEV method and the GeneLab-to-SPOKE pipeline whose
  code this follows (`replication/psev_pipeline.py` is the faithful port).
- Zenodo 10.5281/zenodo.4404618 — precomputed gene PSEVs for the 2019 SPOKE.
- NCBI HomoloGene build 68 — mouse/human orthology.
- GeneLab bulk RNA-seq processing pipeline (GLbulkRNAseq) — source of the RSEM counts and gene annotation.
- Mann HB, Whitney DR (1947); Kolmogorov (1933) / Smirnov (1939) — the two tests, as implemented in SciPy.

---

# Part 2. OSDR-wide run and the baseline-corrected check (added 2026-10-04)

## 13. Which datasets and which samples

`scripts/survey_osdr.py` queries OSDR for every *Mus musculus* study with an RNA-seq assay
(111 studies), caches each study's file listing and ISA metadata, and keeps those with a GeneLab
`RSEM_Unnormalized_Counts` file (106). `scripts/build_manifest.py` then turns the sample tables
into a manifest with one entry per (study, tissue):

- **Control samples.** For spaceflight studies, samples whose `Spaceflight` factor is exactly
  "Ground Control" (basal, vivarium, cohort and rerun controls are not used). For ground-analog
  studies (hindlimb unloading, radiation, diet), samples at the baseline level of every condition
  factor: normally loaded, sham or non-irradiated, not injected, non-treated, normal diet. These
  are labelled `analog control` so they can be filtered out. Genotype is restricted to wild type /
  control. Factors that are not conditions (age, sex, strain, duration, dissection, mission) are pooled.
- **Tissue.** From the `Tissue` factor when a study has several tissues, otherwise from
  `Characteristics[Material Type]`; left/right and regional names are merged (for example "Left
  kidney" → kidney). Multi-tissue studies are split into one entry per tissue (`OSD-457_spleen`).
- **Skipped**, with the reason recorded in `data/_catalog/manifest_skipped.tsv`: cell cultures,
  tumours, plasma, sample-preservation and methods studies with no control group, studies whose ISA
  sample names do not match the counts columns, and entries with fewer than two control samples.

The first 14 datasets (sections 1–12) were run before the manifest existed and keep their original
sample selection (`_GC_` columns).

## 14. Why a baseline correction is needed

With all-positive weights, a sample's score for node *n* is

    score[s, n] = (1/|G|) Σ_g w[g, s] · rank[g, n]

and most of it does not depend on the tissue. Across datasets from different tissues the pooled
Anatomy percentiles correlate at Spearman ≈ 0.96, and 96 % of their variance is one shared profile.
That profile is set by the shape of each node's PSEV distribution across genes, not by expression:
sparsely annotated nodes (hepatic cord, renal papilla, masseter, hippocampal commissure) sit in the
top 1–3 % in every tissue, while heavily annotated organ nodes (liver, kidney, retina) sit near the
bottom in every tissue. The raw check in sections 6–11 therefore mostly measures where a family
happens to fall in this fixed profile, plus a small tissue-specific shift.

## 15. The baseline-corrected score

The tissue-independent part is what the score would be if every gene had the same weight:

    uniform[n] = (1/|G|) Σ_g rank[g, n]

`gc_counts_to_psev.py --compact` computes it in the same pass by carrying one extra pseudo-sample
whose weights are all 1, and stores it with the raw scores in `<name>_GC_ranks.npz`. The corrected
score subtracts it, scaled by the sample's mean weight:

    corrected[s, n] = score[s, n] − mean_g(w[g, s]) · uniform[n]
                    = (1/|G|) Σ_g (w[g, s] − mean_g w[·, s]) · rank[g, n]

which is the covariance, across the dataset's genes, between how strongly a gene is expressed in the
sample and how strongly that gene is linked to node *n*. It is positive when the sample's highly
expressed genes are the ones SPOKE connects to the node, zero when expression and connectivity are
unrelated, and negative when the node's genes are the weakly expressed ones. It uses only the
dataset's own genes and samples; no other dataset enters.

`check_anatomy.py --npz <name>_GC_ranks.npz --corrected` ranks the Anatomy nodes by this corrected
score within each sample (highest = percentile 0) and then runs exactly the same three figures,
tests and verdict as sections 6–11, writing them to `datasets/<name>/baseline_corrected/`.

What the correction does not remove: genes that are highly expressed in every tissue (housekeeping
genes) also tend to be well connected in SPOKE, so nodes linked to ubiquitous genes keep a mild
advantage in every tissue. Removing that needs a reference expression profile, which is a
cross-dataset operation (see the machine-learning notes in the README).

## 16. The compiled tables

`scripts/compile_results.py` reads every `datasets/*/summary.json` (raw) and
`datasets/*/baseline_corrected/summary.json` and writes `results/compiled_verdicts.tsv` (one row per
dataset, raw and corrected columns side by side), `results/tissue_summary.tsv` (one row per tissue)
and `results/tissue_separation_overview*.png`. It tabulates verdicts; it does not test datasets
against each other.

## 17. Ontology-based tissue families

`scripts/ontology_families.py` reads UBERON (`data/_ontology/uberon-basic.obo`, release 2026-10-01)
and defines each family as a root term (liver, kidney, retina, hippocampal formation, telencephalon
for cerebrum, cranial nerve II for optic nerve, ...) plus every term below it through `is_a` and
`part_of`, intersected with SPOKE's Anatomy nodes. `develops_from` is not followed, so primordia are
excluded. `check_anatomy.py --families ontology` uses these masks for the tissue and for the
competition families (minimum 3 nodes; the tissue's own family may be as small as 1) and writes to
`datasets/<name>/corrected_ontology/`. Tests, figures and verdict are unchanged.

## 18. Which Anatomy nodes have gene links, and the linked-node check

The 2019 edge list is not distributed with the PSEVs, so gene links are read from the PSEVs
themselves. `scripts/extract_anatomy_psev.py` stores the 19,567 × 13,257 matrix of gene-PSEV values at
Anatomy nodes. A random walk restarted at gene *g* deposits a floor of about 0.09 / degree(*g*) on each
direct neighbour, so a node linked to a gene has a visibly higher value than one reached only through
other nodes. 74 Anatomy nodes exceed 1e-4 for more than 99 % of genes; the 75th (gall bladder) does so
for 6 %. Those 74 are the linked set. A node's mean PSEV value across genes is used as its link
weight (in an undirected walk the stationary mass is proportional to degree).

`scripts/check_linked.py` ranks the 74 linked nodes by baseline-corrected score within each sample
(1 = top), averages the ranks over samples, and reports the rank of the tissue's own linked node(s):
the linked nodes inside its ontology family, or, if there are none, the nearest linked ancestors of
the root term (flagged proxy). CLEAR = a target ranked 1st; PARTIAL = best target in the top 7;
NONE otherwise; NO LINKED NODE when neither the family nor any ancestor is linked. Figures:
`A_linked_nodes` (top 20 of the 74 with the targets highlighted) and `C_per_sample`.

## 19. Stage-1 pilot: four ways of scoring (2026-10-11)

Ten datasets, one per tissue (retina, kidney, heart, spleen, liver, adrenal gland, lung,
cerebellum, white adipose tissue, blood), scored four ways and judged on the gene-linked nodes
(section 18) and on ontology families (section 17). All arms are baseline-corrected and each dataset
is processed alone. `scripts/pilot_stage1.py` runs the checks and tabulates `results/pilot_stage1.tsv`.

| arm | what changes | flag |
|---|---|---|
| A | nothing: within-gene ranks, log2(CPM+1) weights (the current method) | — |
| Z | the gene's z-scores are weighted directly instead of their within-gene ranks | `--zscore-arm`, `--score-key zscore` |
| B | observed corrected score converted to a z against 50 permutations of the sample's weights among genes with similar hub rank (20 quantile bins of `hub_rank` from `gene_hub_rank.py`) | `--perm-cov`, `--score-key perm` |
| D | housekeeping genes removed before weighting (HRT Atlas mouse list, 3,174 Ensembl ids; about a quarter of expressed genes) | `exclude_genes` in the manifest |

`hub_rank` is, for each gene, the mean within-gene rank it assigns to the 74 linked Anatomy
nodes, computed once over all 19,567 genes; `n_high_z` (nodes with z > 3) is stored alongside.
Permuting weights within hub-rank bins keeps each sample's overall tendency to push the hub nodes
up while breaking the gene-specific pairing, so a node's permutation z reflects only the
tissue-specific part.

Success measures fixed before running: number of tissues whose own linked node ranks first (and in
the top three) of 74; Spearman between a linked node's link weight and its rank in its own tissue
(−0.85 under arm A; the goal is near zero); ontology-family verdicts.

**Result.** Own-node ranks (of 74) were unchanged within a few places across all four arms for the
failing tissues (liver 48–57, lung 60–69, cerebellum 64–68, adipose 63–65, adrenal cortex 34–42) and
the link-weight relationship stayed at Spearman −0.87 to −0.93. Arm D equals arm A; arm Z lifts
blood plasma (13 → 3) and some family medians but not the organ nodes; arm B does not help. Arm A
is retained. The pilot projections (`psev_out/*__pilot_GC_ranks.npz`, with `scores_z`, `uniform_z`
and `perm_scores`; `*__nohk_*`) and checks (`datasets/*__pilot/`, `datasets/*__nohk/`) are kept.

## 20. Stage 2: other SPOKE builds

The same scripts run on any gene-PSEV build directory: `gc_counts_to_psev.py --zip <dir> --spoke <dir>`
(the build's own `node_info.tsv` defines the node universe; 2025 builds have 547,873 nodes),
`extract_anatomy_psev.py --build <dir> --spoke <dir> --out psev_out/_anatomy_psev_<tag>.npz` to find
that build's gene-linked Anatomy nodes, `check_linked.py --spoke <dir> --anatomy-psev ... --tag <tag>`
and `check_anatomy.py --corrected --families ontology --spoke <dir>`. `pilot_stage2.py --tag <tag>
--build <dir> --run-checks` runs both for the pilot datasets and writes `results/pilot_stage2_<tag>.tsv`;
`plot_stage2.py` draws the build comparison. The 2025 Anatomy vocabulary with gene links is GTEx-like
(right lobe of liver, upper lobe of left lung, skin of leg), so the ontology families pick them up
through `part_of` without any change to the family definitions.

## 21. Full run on the 2025 'upregulates-only' build

`run_upauG_all.sh` projects all 104 control sets (`data/_catalog/manifest_2025upAuG.json`: the 90 manifest
entries plus the first 14 with their original `_GC_` column selection) on
`/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only` in three batches and runs
`check_build.py` (linked-node check with `psev_out/_anatomy_psev_2025upAuG.npz`, 44 linked nodes, and the
ontology-family check against that build's node table) after each. `compile_build.py --tag 2025upAuG`
tabulates them. Outputs sit beside the 2019 ones as `datasets/<name>__2025upAuG/`.

## 22. Cell-type route

`scripts/celltype_check.py --tag <tag> --build <dir>` takes the baseline-corrected scores already stored
in `<name>_GC_ranks.npz` (all 547,873 nodes), keeps the CellType columns, ranks them per sample (1 = top),
pools the ranks over samples, and reports the best-ranked cell type whose name matches a tissue regex
(`FAM` in the script) together with the top five cell types. No projection is re-run. Cell types reach
the walk through their own Gene edges (59 cell types) and through `PARTOF_CTpA` links to Anatomy nodes.

## 23. Linked set from the edge list

When the live graph is reachable, `inspect_aug_coverage.py` writes the per-Anatomy counts of the three
Anatomy–Gene edge types. `check_linked.py --linked-ids FILE` then uses an explicit UBERON id list as the
linked set instead of the PSEV-threshold detection of section 18, which under-counts nodes with few
edges. The 60-node set (≥ 200 `UPREGULATES_AuG` edges) is the default for the upregulates-only build;
the 72-node set (≥ 1 edge) is kept for reference (`corrected_linked_2025upAuG72/`).

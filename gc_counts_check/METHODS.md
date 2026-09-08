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

# Presenter guide: results/SPOKE_tissue_validation.pptx

One-sentence story: **when we project a tissue's control-sample gene counts onto SPOKE, does the tissue's own
Anatomy node come out on top? With the full graph, no. Once we found out why and rebuilt the PSEVs with only the
up-regulation Anatomy–Gene edges, yes for almost every tissue we can test.**

Each slide: what to say, what to point at, and the question you are most likely to get.

---

## Slide 1 – Title
**Say:** "This talk is about a sanity check. Before we use SPOKE embeddings of OSDR data to study spaceflight, we
need to show that the embedding knows what tissue a sample came from. I ran that check on every mouse control
dataset on OSDR: 104 control sets, 31 tissues, 916 samples, on the current 2025 SPOKE."

## Slide 2 – Why control samples, and why anatomy
**Say:** "Control animals never flew, so the only thing their samples should 'say' when projected onto SPOKE is
their tissue. That makes them a built-in positive control. If a liver control doesn't point at the liver node, we
can't trust what it says about spaceflight."
**Point at:** the tiles. SPOKE has 13,895 Anatomy nodes; only 60 have enough gene edges to be tested directly.
**Stress:** every dataset is judged alone. We never compare one dataset against another, because batch effects
between studies would swamp everything.
**Likely question:** "What's a control set?" One study's ground-control (or analog baseline) samples for one
tissue. A study with three tissues gives three control sets.

## Slide 3 – Pipeline
Walk the five boxes, one sentence each:
1. **Counts → weights.** "Normalize each sample to counts per million, drop genes with CPM below 1, weight each
   gene by log2(CPM+1). Highly expressed genes count more, but not thousands of times more."
2. **Mouse → SPOKE genes.** "SPOKE genes are human, so mouse Ensembl ids go to human Entrez ids through HomoloGene."
3. **PSEV projection.** "Each gene has a PSEV: a vector over all 548k SPOKE nodes saying how close each node is to
   that gene by random walk. For each gene we rank the nodes, then take the weighted average rank across the
   sample's genes. Nodes with a low average rank are the ones the sample's expression 'points at'."
4. **Baseline correction.** "The step we added; slide 5 explains why."
5. **Rank & test.** "Three questions, all inside the dataset: does the tissue's own node rank first among the 60
   gene-linked Anatomy nodes; does the tissue's whole UBERON family sit near the top (QQ plot, Mann-Whitney, KS);
   and for tissues with no linked node, does the tissue's cell type rank first."
**Likely question:** "Why ranks and not raw PSEV values?" PSEV values span orders of magnitude and are dominated
by hub nodes; ranks make genes comparable.

## Slide 4 – Section divider
**Say:** "Three things we learned along the way changed the configuration."

## Slide 5 – Finding 1: the raw ranking is mostly the graph
**Say:** "When I first ran this, every tissue gave almost the same Anatomy ranking. A liver dataset's raw ranking
correlates at 0.93 with the ranking you get when every gene is weighted equally, and datasets from seven
different tissues agree with each other at 0.87. The same five kidney-cortex nodes head every one of them. That's
not biology; it's the graph. Because every weight is positive, a node's score is mostly how central it is."
**Fix:** "So within each dataset we also score a pretend sample where every gene has weight 1. That's the pure
graph baseline. We subtract it. What's left is only the part of the score that depends on which genes are high
in this sample."
**Point at the table:** "OSD-245 liver: right lobe of liver goes from 161st to 9th of all 13,895 Anatomy nodes,
and from 4th to 1st among the 60 linked nodes. Cortex of kidney drops from 1st to 2,882nd."
**Point at the QQ plot:** "The orange curve is the liver family of 70 UBERON nodes. The x-axis is where they'd
fall if they were random; the y-axis is where they actually fall, 0 is the top. Hugging the bottom means most of
the liver family is in the top percent of Anatomy nodes. The grey band is what random node sets do; grey lines
are other organ families in the same dataset. Thin lines are individual samples, all agreeing."
**Likely question:** "Isn't subtracting a baseline just fitting to the result?" No: the baseline uses the graph
alone, no expression data, and the same correction is applied identically to every dataset.

## Slide 6 – Finding 2: few Anatomy nodes carry gene links
**Say:** "SPOKE has 13,895 Anatomy nodes, but only 92 have any Anatomy–Gene edge, 72 have up-regulation edges, and
60 have enough of them to test. Those edges come from human atlases that sampled a few dozen tissues. Every other
node, including thymus, retina and bone marrow, only gets a score second-hand, through ontology edges."
**Consequence:** "So the honest test isn't 'is liver in the top 1 % of 13,895', it's 'among the nodes that can
actually be tested, where does liver rank'. On the full graph the answer was 35th of 92, and the
pattern was systematic: a node's rank in its own tissue was predicted by how many edges it had, correlation
−0.95. The lightly linked nodes won everywhere."
**Likely question:** "Why would the number of edges decide the rank?" Preview slide 7: the walk treats a blanket
'expresses' edge the same as a specific 'upregulates' edge.

## Slide 7 – Finding 3: edge type is the cause
**Say:** "PSEVs are built on a binary adjacency: an edge is an edge. But the graph has 1.58 million 'expresses'
edges, where for example blood is linked to 18,879 of 19,072 genes. That's not information, it's a blanket. Only
the 86,000 'upregulates' edges, which require log2 fold-change of at least 1, say something specific."
**Point at the figure** (each row a tissue, x = how far down its own node ranks, one marker per edge set):
- Squares, all edges: liver, spleen, blood far right.
- Triangles, no 'expresses' edges: worse everywhere.
- Inverted triangles, housekeeping genes removed on the full graph: same as all edges. Not a gene problem.
- **Diamonds, 'upregulates' only: every tissue that has a node is at the left edge, top 3.** The connectivity
  correlation drops to about zero.
- Open markers at the far right mean that edge set has no linked node for the tissue. Retina, kidney, heart lose
  their node when only up-regulation edges are kept, which is why slide 11 has a cell-type route.
**Say:** "Same samples, same scoring. The only thing that changed is which edges the random walk could use."

## Slide 8 – Section divider
**Say:** "So the validated configuration is the same pipeline on a different edge set."

## Slide 9 – Validated configuration table
Read it row by row. Key points: the graph is the 2025-03 SPOKE with PSEVs rebuilt using only the 85,530
UPREGULATES Anatomy–Gene edges, everything else in SPOKE unchanged. CLEAR means the tissue's own node ranks 1st of
60, PARTIAL means top 7. The cell-type route covers tissues without a linked Anatomy node. The last row lists the
three scoring tricks that were tested and rejected, including housekeeping-gene removal.

## Slide 10 – Results on 104 datasets
**Point at the table:** 50 CLEAR, 32 PARTIAL, 2 NONE, 20 with no linked node.
**Point at the figure:** "Each row is a tissue, number in parentheses is how many datasets. Each marker is one
dataset, placed at the rank of its own node. Green at 1 is a win. Liver 14 of 14, all five hindlimb muscles 20
of 20, cerebrum, brain, white adipose and blood every time; kidney 4 of 6 with the other two second."
**Explain the near-misses:** "The yellow triangles are mostly losing to a biologically related node: spleen to
blood in 6 of 8, heart to gastrocnemius in all 6, adrenal gland and skin to the aorta, lung to arteries. Rank 2 to
6 of 60, not failures."
**Red crosses:** one colon set and brown adipose. The open squares at the far right are the 20 datasets whose
tissue has no node among the 60; those go to the cell-type route.
**Also say:** ground-analog baselines (hindlimb unloading, radiation controls) behave like spaceflight ground
controls.

## Slide 11 – Tissue coverage in three groups
**Say:** "20 tissues validate directly through their Anatomy node. Four more validate through a cell type: retina
through photoreceptor cell, rank 1 in four of five datasets and 2 in the fifth; bone marrow through monocyte;
placenta through syncytiotrophoblast; cecum through Paneth cell. Four are not testable at all: bone, mandible,
temporal bone, optic nerve, because no node with gene edges exists for them."
**Bottom line:** "27 of 31 tissues have a direct within-dataset test. Thymus, eye, colon and brown adipose pass
only weakly and are labelled that way."
**Likely question:** "Why no thymus?" The gene edges come from human atlases like GTEx, which never sampled
thymus, retina, bone marrow or placenta. A data-coverage gap in SPOKE, not a method failure.

## Slide 12 – Liver example
**Say:** "One concrete example. OSD-245, 19 liver ground controls, same 2025 SPOKE on both sides. With all
Anatomy–Gene edges the liver node ranks 35th of 92 linked nodes. With only the up-regulation
edges, right lobe of liver is 1st of 60, followed by adrenal gland and kidney: metabolic organs, which is
sensible. Same samples, same weights, same scoring; only the edge set changed."

## Slide 13 – The same ten tissues on four edge sets
**Say:** "Here is the comparison as a table: one row per pilot tissue, one column per edge set, cells are the
rank of the tissue's own node. The bold column is the up-regulation-only set. The bottom row is the correlation
between a node's edge count and its own-tissue rank: strongly negative with all edges, about zero with
up-regulation edges only. Removing housekeeping genes on the full graph changes nothing, which is why we stopped
looking at the gene weights and rebuilt the graph."
**Why this slide matters:** it pre-empts "did you try X?" and justifies the rebuild, which was the expensive step.

## Slide 14 – Limitations
One line each:
- Coverage follows the human atlases: 60 Anatomy nodes, 59 cell types. Mouse tissues are validated via human
  counterparts.
- Residual hubs: gastrocnemius is in the top four of 43 datasets; sinoatrial myocytes recur in cell-type rankings.
- UBERON families still mix organ tissue with its vessels and embryonic parts for some organs.
- Small sets: some tissues have one dataset; OSD-270 heart has two samples; OSD-289 GC includes vivarium controls.
- RR-1 thymus and liver counts share a technical factor, so those two datasets are not independent.
- Direction is lost: 'downregulates' edges had to be dropped because the walk is unsigned.

## Slide 15 – Next steps
1. **Freeze:** configuration and tissue table are reported; the standalone repository runs the check on any
   counts table.
2. **Benchmark:** a tissue classifier on the corrected SPOKE vectors versus one on raw expression, held out by
   study. It tests two things the anatomy check cannot: transfer across studies, and whether tissue identity
   survives in the Pathway and Disease parts of the vector that the flight analysis will read.
3. **Flight:** project flight samples with the identical pipeline, confirm the tissue still ranks first (a quality
   gate), then test flight vs ground within study and tissue, checking known effects before new ones.
4. **Graph:** background experiment with a stricter edge threshold (log2 fold-change ≥ 2); mouse atlas edges only
   if thymus or bone become important.

---

## Three answers to have ready

**"What is a PSEV in one sentence?"** A vector over every SPOKE node giving how often a random walker that keeps
restarting at one gene lands on each node; it encodes the gene's neighbourhood in the whole graph.

**"What exactly does CLEAR mean?"** After baseline correction, the tissue's own Anatomy node has the best mean
rank, across that dataset's control samples, of the 60 Anatomy nodes that have at least 200 up-regulation edges.

**"Why does this matter for spaceflight?"** Any flight-vs-ground difference in SPOKE space is only interpretable
if the dominant signal in SPOKE space is biology. Showing the tissue signal is recovered within each dataset is
the minimum evidence that it is.

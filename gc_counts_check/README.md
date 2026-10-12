# GC Counts Check — does a ground-control expression profile point at its own tissue in SPOKE?

Each OSDR dataset is judged **on its own**: its ground-control (GC) counts are projected onto the
2019 SPOKE gene PSEVs (Nelson et al., Life 2021; Zenodo 10.5281/zenodo.4404618), and the resulting
Anatomy-node ranking is checked for separation of the known source tissue. Nothing is compared
across datasets.

## Folder layout

```
scripts/gc_counts_to_psev.py   counts -> per-sample SPOKE rank table (2019 PSEVs), any dataset via --counts
scripts/check_anatomy.py       single-dataset anatomy check -> datasets/<name>/ (3 figures + stats + verdict)
scripts/families.py            Anatomy node families (UBERON name regexes) used by the check
scripts/style.py               figure style (white background, fixed tissue colours)
run_all.sh                     projection + one check per dataset for the first 14 datasets
run_osdr_wide.sh               OSDR-wide run (survey -> manifest -> projection -> checks -> compile)
data/                          GeneLab RSEM unnormalized counts + ISA metadata (OSD-564, 244, 245, 246, 288, 289)
psev_out/                      rank tables (Rank_by_type_<sample> columns), prepared gene weights, run log
datasets/<OSD>/                A_qq_plot, B_organ_families, C_per_sample (.png/.svg), summary.json,
                               organ_families.tsv, per_sample.tsv
archive/                       earlier cross-dataset analysis (see below); not part of the check
```

## Datasets (GC samples only, `_GC_` columns)

| dataset | tissue | GC samples | note |
|---|---|---|---|
| OSD-564 | hippocampus | 5 | RR-10 |
| OSD-244 | thymus | 19 | RR-1; 9 LAR + 10 ISS-T (carcass-dissected) |
| OSD-289 | thymus | 6 | MHU-1/2; the 3 MHU-2 "GC" columns are vivarium controls per ISA |
| OSD-246 | spleen | 16 | RR-1 |
| OSD-288 | spleen | 3 | MHU-1 |
| OSD-245 | liver | 19 | RR-1 |

OSD-4 (Affymetrix thymus, no counts) was excluded at the user's request.

RR-3 (SpaceX-8) payload, added 2026-09-08 (`gc_counts_to_psev.py --only RR3`):

| dataset | tissue | GC samples | note |
|---|---|---|---|
| OSD-137 | liver | 6 | |
| OSD-161 | adrenal gland | 5 | |
| OSD-162 | eye | 5 | whole eye |
| OSD-163 | kidney | 6 | left kidney |
| OSD-168 | liver | 4 | RR-1 + RR-3 study; only the `_RR3_GC_` (with ERCC spike-in) columns used |
| OSD-194 | retina | 4 | left/right retina |
| OSD-270 | heart | 2 | |
| OSD-352 | brain | 3 | bulk RNA-seq columns of a multi-modal study |

## Running the check on any dataset

```bash
# from an existing rank table
python3 scripts/check_anatomy.py --ranks psev_out/OSD-246_GC_counts_ranks_and_rank_by_type_for_meta0_1.tsv --tissue spleen --name OSD-246
# from a raw GeneLab counts file (runs the ~7 min projection first; needs psev_data/gene_psev.zip)
python3 scripts/check_anatomy.py --counts path/to/GLDS-xxx_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv --tissue liver --name OSD-xxx --samples _GC_
```

`--tissue` is a built-in family (hippocampus, thymus, spleen, liver, heart, kidney, lung, skeletal
muscle, bone marrow, skin, intestine, testis, blood, cerebral cortex) or use `--family-regex` /
`--exclude-regex` for any other UBERON term. `--out` changes the destination (default `datasets/`).

## Method

A step-by-step guide with the statistics explained is in `METHODS.md`.

1. **Projection.** Per gene and GC sample the weight is log2(CPM+1) (genes with CPM ≥ 1 in ≥ 1 sample);
   ENSEMBL → mouse Entrez (GeneLab annotation) → human Entrez (HomoloGene build 68) → SPOKE v2 gene
   node. Same maths as `replication/psev_pipeline.py` (z-score each gene PSEV over the dataset's genes,
   rank nodes within the gene, expression-weighted mean, rank; `Rank_by_type_` per sample).
2. **Anatomy check** on the 13,257 Anatomy nodes, percentile of rank-by-type (0 = most associated with
   the sample's expressed genes), pooled over the dataset's samples:
   - **A. QQ plot** — sorted percentiles of the tissue family's nodes vs the percentiles a random set of
     Anatomy nodes would have; 95 % envelope from 2,000 random sets of the same size; one-sided KS test
     and Mann-Whitney U (family vs all other Anatomy nodes).
   - **B. Organ families** — the family's median percentile against the other organ families in the
     same dataset (black bar = median); rank among families.
   - **C. Per sample** — each sample's family median and its own Mann-Whitney p.
3. **Verdict** (thresholds are flags): CLEAR = pooled MWU p < 0.01, KS p < 0.01, ≥ 80 % of samples
   individually p < 0.05, family ranked first; PARTIAL = pooled p < 0.05 but a criterion fails;
   NONE = no pooled shift.

## Results (2026-09-08)

| dataset | tissue | verdict | pooled MWU p | KS p | samples p<0.05 | family rank |
|---|---|---|---|---|---|---|
| OSD-564 | hippocampus | CLEAR | 1.8e-18 | 3e-24 | 5/5 | 1/13 |
| OSD-246 | spleen | CLEAR | 1.6e-3 | 2e-5 | 16/16 | 1/13 |
| OSD-244 | thymus | PARTIAL | 0.024 | 3.3e-5 | 18/19 | 2/13 (behind the 5-node bone marrow family) |
| OSD-288 | spleen | PARTIAL | 4.0e-3 | 1e-4 | 3/3 | 4/13 |
| OSD-289 | thymus | NONE | 0.11 | 4e-3 | 0/6 | 5/13 |
| OSD-245 | liver | NONE | 0.44 | 0.13 | 2/19 | 8/13 |

RR-3 payload (families for adrenal gland, eye, retina and brain added to `families.py`; the
competition set stays the fixed 13 organ families plus the dataset's own family):

| dataset | tissue | verdict | pooled MWU p | KS p | samples p<0.05 | family rank | family nodes in top 2.5 % (expected) |
|---|---|---|---|---|---|---|---|
| OSD-163 | kidney (167 nodes) | PARTIAL | 1.3e-7 | 1.1e-8 | 6/6 | 7/13 | 7 (4.2) |
| OSD-352 | brain (578 nodes) | PARTIAL | 3e-78 | 2e-93 | 3/3 | 3/14 (behind hippocampus, cerebral cortex) | 17 (14.5) |
| OSD-270 | heart (97 nodes) | PARTIAL | 0.037 | 0.056 | 2/2 | 7/13 | 1 (2.4) |
| OSD-137 | liver (85 nodes) | NONE | 0.66 | 0.36 | 0/6 | 8/13 | 1 (2.1) |
| OSD-168 | liver (85 nodes) | NONE | 0.54 | 0.38 | 0/4 | 8/13 | 0 (2.1) |
| OSD-161 | adrenal gland (31 nodes) | NONE | 0.32 | 0.0044 | 0/5 | 12/14 | 2 (0.8) |
| OSD-162 | eye (169 nodes) | NONE | 1 | 0.8 | 0/5 | 14/14 | 1 (4.2) |
| OSD-194 | retina (37 nodes) | NONE | 1 | 0.87 | 0/4 | 14/14 | 0 (0.9) |

Reading notes:
- RR-3 shows the same pattern as the first set: the tissue's parenchymal core sits at the very top
  (liver lobule 2.6 %, zona fasciculata 1.0 % and zona glomerulosa 1.5 % of the adrenal, corneal
  epithelium 2.3 % and ciliary body 2.6 % of the eye, renal papilla 0.9 %, cardiac Purkinje fiber 1.5 %,
  olfactory bulb layers 0.6 %), while the broad name-based family is diluted by vessels, ducts and
  primordia, so the family-level verdict is sensitive to how the family is drawn.
- Retina is the one real failure: all 37 retina nodes sit at the bottom of the ranking (median 93 %)
  in every retina sample, i.e. retina anatomy is anti-associated with retina expression on the 2019
  build. The eye dataset shows the same for its retina members. This looks like a property of how
  retina is wired in SPOKE v2 rather than of the samples.
- OSD-270 has only two GC samples; its verdict rests on very little.
- Thymus and spleen families are bimodal: organ/cortex/medulla/pulp terms sit near the 10th percentile,
  vessels and embryonic primordia at the bottom. The families separate on their parenchymal core.
- The 85-node liver family mixes parenchyma with bile-duct, ligament and vessel terms; the parenchymal
  nodes are at the top (hepatic cord 2.5 %, liver lobule 2.6 %, hepatic acinus 3.0 %) while the rest sit
  at the bottom, so the broad family reads as NONE. `--tissue "liver parenchyma"` (24 nodes, post hoc)
  gives median 28 %, MWU p 8e-6.
- On the 2019 build the absolute top of every Anatomy ranking is UBERON life-stage pseudo-nodes and
  brain/neural hub nodes regardless of tissue (listed under `top_real_anatomy` in each summary.json), so
  the check is about where the family sits, not what is at rank 1.
- Expression-side QC (archived) confirmed every dataset expresses its own tissue signature, so failures
  here are properties of the projection, not of the counts. It also found a strong technical factor
  shared by the RR-1 thymus and liver replicates (samples pair by replicate index, residual r ≈ 0.99).

## OSDR-wide run (2026-10-04): 104 datasets, 31 tissues, 916 control samples

Every mouse bulk RNA-seq study on OSDR with GeneLab counts and a usable control group was added
(`run_osdr_wide.sh`; about 4.5 h of compute for 104 datasets on an 8-core laptop):

```
scripts/survey_osdr.py      111 mouse RNA-seq studies -> data/_catalog/survey.tsv (+ ISA metadata, file listings)
scripts/build_manifest.py   control samples per (study, tissue) -> data/_catalog/manifest.json, manifest_skipped.tsv
scripts/gc_counts_to_psev.py --manifest ... --compact
                            -> psev_out/<name>_GC_anatomy_ranks.tsv (raw Anatomy ranks)
                               psev_out/<name>_GC_ranks.npz (ranks, raw scores, uniform-weight profile, mean weights; all 389,297 nodes)
scripts/check_manifest.py [--corrected]     one check per dataset -> datasets/<name>/ and datasets/<name>/baseline_corrected/
scripts/compile_results.py                  -> results/compiled_verdicts.tsv, tissue_summary.tsv, tissue_separation_overview*.png
```

83 entries are true ground-control groups from spaceflight missions and 21 are baseline controls of
ground-analog studies (`control_type` column). Selection rules are in METHODS.md section 13.

### The raw check mostly measures a tissue-independent baseline

Across datasets of different tissues the pooled Anatomy percentiles correlate at Spearman 0.96;
96 % of their variance is one shared profile. Nodes such as hepatic cord, renal papilla, masseter
and hippocampal commissure sit in the top 1–3 % in every tissue, and the organ nodes liver, kidney
and retina sit near the bottom in every tissue. Two earlier statements in this README are therefore
wrong as written: the "parenchymal core at the top" pattern is the baseline, not tissue signal, and
retina is not anti-associated. `METHODS.md` sections 14–15 describe the within-dataset correction
(subtract the uniform-weight profile), which every dataset now has in `baseline_corrected/`.

### Results per tissue (baseline-corrected; raw in `results/tissue_summary.tsv`)

| verdicts | raw | baseline-corrected |
|---|---|---|
| CLEAR | 6 | 20 |
| PARTIAL | 33 | 64 |
| NONE | 65 | 20 |
| pooled Mann-Whitney p < 0.01 | 35 | 79 |
| tissue family ranked first of the organ families | 6 | 20 |

- **Separate clearly** (family median in the top ~15 %, usually ranked first): retina (5/5 CLEAR;
  the `retina` node is ranked 1 of 13,257), bone marrow, cecum, kidney, heart, hippocampus, brain,
  colon, cerebrum, mammary gland.
- **Separate, but not ranked first** (PARTIAL in every dataset): spleen (8), thymus (5), all five
  skeletal muscles (20 datasets; top nodes are cardiovascular and striated-muscle terms), skin (7),
  eye, placenta. Liver is borderline: 10 PARTIAL, 4 NONE across 14 datasets, family median ≈ 31 %,
  and the top nodes of a liver sample are kidney and renal terms.
- **Do not separate**: adrenal gland (4/4), cerebellum (3/3; brain samples all point at forebrain),
  lung (2 of 3), white and brown adipose, bone tissue, mandible, temporal bone, optic nerve, whole
  blood. These are candidates for weak or missing tissue wiring in SPOKE v2.
- Analog-study controls behave like ground controls for the same tissue.
- Residual artifact: a few nodes (for example `mammary gland smooth muscle`, `cardiovascular system`)
  rank near the top in many unrelated tissues after correction, because they are linked to genes
  that are highly expressed everywhere. A reference-centred weight would remove this.

### Using the data for machine learning

Per-sample vectors for all 389,297 nodes are in `psev_out/<name>_GC_ranks.npz`
(`scores - mean_weight[:, None] * uniform` is the corrected embedding). Use the corrected vectors,
split by study, compare against a gene-expression baseline, and reduce dimensions (one node type, or
principal components) before fitting.

## Two no-re-projection refinements (2026-10-05)

Both run on the stored baseline-corrected scores (`scripts/check_manifest.py --variant ontology|linked`).

**1. Ontology families** (`scripts/ontology_families.py`, `check_anatomy.py --families ontology`,
output `datasets/<name>/corrected_ontology/`). A family is a UBERON root term plus everything below
it by is_a / part_of, instead of a name regex. Result: sharper, not more permissive. 17 CLEAR /
58 PARTIAL / 29 NONE (name families: 20 / 64 / 20). Kidney becomes CLEAR in 6 of 6 and mammary
gland CLEAR; retina's family median drops to 0.2 %. But the true organ families of liver, adrenal
gland, lung and cerebellum sit at the *bottom* of the ranking (medians 73 %, 99.7 %, 92 %, 78 %); the
name-based liver family had looked borderline only because it included non-liver terms.

**2. Gene-linked nodes only** (`scripts/extract_anatomy_psev.py`, `scripts/check_linked.py`, output
`datasets/<name>/corrected_linked/`). Reading the Anatomy columns of all 19,567 gene PSEVs shows that
only **74 of the 13,257 Anatomy nodes are linked to genes** in the 2019 graph (a linked node has a
PSEV value above 1e-4 for essentially every gene; the 75th node has it for 6 %). Spleen, thymus,
skeletal muscle, colon, cecum, hippocampus, cerebral cortex and all bone terms have no gene links;
their signal is second-hand through neighbours such as lymphoid tissue, hematopoietic system and
digestive system. The check ranks the 74 linked nodes per sample and reports where the tissue's own
node falls (or its nearest linked ancestor, flagged "proxy").

| tissue | linked node(s) | rank of 74 in own samples |
|---|---|---|
| retina | retina | 1–2 |
| brain | forebrain, telencephalon | 1–2 |
| kidney | kidney, nephron tubule | 2–3 |
| heart | myocardium | 2–3 (the `heart` node itself: 65) |
| spleen, thymus | hematopoietic system (proxy) | 2 |
| colon, cecum | digestive system (proxy) | 1–2 |
| skin | skin of body | 11 |
| mammary gland, blood | own node | 12, 13 |
| bone marrow, placenta | bone marrow, decidua | 26–29 |
| adrenal gland | adrenal cortex / adrenal gland | 34 / 72 |
| liver, adipose, cerebellum, lung | own node | 60–70 |
| skeletal muscles, bone | none, and no linked ancestor | not testable |

**Why the failures fail.** Across the 31 linked nodes that are some dataset's own tissue, how heavily a
node is linked to genes predicts how badly it ranks in its own samples: Spearman −0.85 (p 1.5e-9).
The most heavily linked nodes (brain, adrenal gland, heart, lung, cerebellum, liver) are connected to
so many genes that the link carries no tissue information: for liver, the correlation across genes
between expression in liver samples and the gene's PSEV value at the liver node is 0.00, against
+0.33 for retina and +0.36 for kidney. Lightly linked nodes (retina, kidney, myocardium, forebrain)
are selective and work. The next improvement should therefore correct for connectivity (a
degree-matched permutation null, or normalising each node's PSEV by its link weight); removing
housekeeping genes alone is unlikely to rescue liver, lung, adrenal or cerebellum.

## Stage-1 pilot: can reweighting rescue the failing tissues? (2026-10-11)

Ten tissues, four scoring arms (METHODS.md section 19), judged by where the tissue's own gene-linked
node ranks among the 74 linked nodes. `results/pilot_stage1.tsv`.

| tissue (own node) | A rank-based | Z z-score | B permutation z | D housekeeping removed |
|---|---|---|---|---|
| retina | 1 | 1 | 1 | 1 |
| kidney | 1 | 2 | 2 | 1 |
| heart (myocardium) | 2 | 4 | 15 | 3 |
| blood (blood plasma) | 13 | 3 | 5 | 11 |
| adrenal gland (adrenal cortex) | 37 | 36 | 42 | 34 |
| liver | 57 | 52 | 48 | 52 |
| white adipose tissue | 65 | 63 | 63 | 64 |
| cerebellum | 68 | 66 | 64 | 67 |
| lung | 69 | 69 | 60 | 67 |
| Spearman(link weight, own rank) | −0.88 | −0.93 | −0.87 | −0.87 |

No arm moves liver, lung, adrenal gland, cerebellum or adipose out of the bottom third, and the
connectivity relationship is unchanged in every arm. Housekeeping removal (D) is indistinguishable
from the current method; the z-score projection (Z) helps blood and softens some family medians
(liver family 57 → 19 %, cerebellum 80 → 27 %) without moving the organ nodes themselves; the
permutation null (B) does not help and costs heart. Decision: keep the rank-based, baseline-corrected
score (A) as the method. The failures are a property of the 2019 graph's edges for those nodes, not
of the weighting, so Stage 2 tests the graph: the 2025 SPOKE build for the failing tissues, and a
no-propagation control on the 74 linked nodes.

## Stage 2: the same ten tissues on the 2025 SPOKE builds (2026-10-11)

Three gene-PSEV builds of the 2025-03 SPOKE snapshot exist on the external drive
(`/Volumes/justinytlin/SPOKE_psev/`, see `psev_build/README.md`): the full graph; the graph without the
1.58 M `EXPRESSES_AeG` edges (keeping 494 k `DOWNREGULATES_AdG` + 86 k `UPREGULATES_AuG`); and the graph
with only the 86 k `UPREGULATES_AuG` Anatomy–Gene edges. The pilot datasets were projected on each
(`pilot_stage2_<tag>.json`, `gc_counts_to_psev.py --zip <build> --spoke <build>`), the gene-linked
Anatomy set was read from each build (`extract_anatomy_psev.py --build`), and the linked-node and
ontology checks were run (`pilot_stage2.py`). `results/pilot_stage2_<tag>.tsv`,
`results/pilot_stage2_builds.png`.

Rank of the tissue's own gene-linked node (number of linked nodes in brackets):

| tissue | 2019 (74) | 2025 all edges (92) | 2025 no 'expresses' (73) | 2025 'upregulates' only (44) |
|---|---|---|---|---|
| retina | 1 | no node | no node | no node |
| kidney | 1 | 21 | 38 | no node |
| heart | 2 | 7 | 45 | no node |
| spleen | proxy 2 | 75 | 48 | 2 |
| liver | 57 | 35 | 45 | **1** |
| adrenal gland | 37 | 5 | 38 | **2** |
| lung | 69 | 11 | 30 | **3** |
| cerebellum | 68 | 1 | 20 | **1** |
| white adipose | 65 | 10 | 37 | **1** |
| blood | 13 | 82 | 73 | **1** |
| Spearman(link weight, own rank) | −0.88 | −0.91 | −0.59 | −0.04 |

Housekeeping-gene removal on the full 2025 graph (`2025nohk`) reproduces the full-graph column within
a few places (kidney 20, liver 31, spleen 74, blood 79, same universal-top nodes), as it did on 2019.
| nodes in the top 3 of ≥ 3 datasets | smooth muscle tissue, cardiovascular system | skeletal muscle tissue, cerebellum, corpus callosum, popliteal artery | putamen, amygdala, blood | upper lobe of left lung only |

- **The full 2025 graph has the same defect as 2019, inverted.** Link weight still predicts own-tissue
  rank (−0.91), but now the most lightly linked nodes (cerebellum, skeletal muscle tissue, corpus
  callosum, popliteal artery) top eight or nine of the ten datasets whatever the tissue, and the most
  heavily linked (blood, heart left ventricle) sink. Cerebellum's "CLEAR" there is largely that artifact.
- **Removing the 'expresses' edges but keeping 'downregulates' fails everywhere**: the top nodes become
  putamen and amygdala in every dataset.
- **Keeping only 'upregulates' edges removes the artifact.** Every testable tissue's own node is in the
  top three, the top nodes differ by tissue (right lobe of liver, left adrenal gland, cerebellar
  hemisphere, subcutaneous adipose tissue, blood), the connectivity correlation is gone (−0.04) and no
  node tops many datasets. The cost is coverage: only 44 Anatomy nodes carry up-regulation edges, so
  retina, kidney and heart have no testable node on that build (their ontology families still read
  PARTIAL / PARTIAL / CLEAR through neighbours).

**Interpretation.** The PSEV random walk uses a binary adjacency, so an Anatomy–Gene edge counts the
same whether it says the gene is expressed in, up-regulated in, or down-regulated in the tissue. The
blanket "expresses" edges (2019: nearly every gene for 74 tissues; 2025: 1.58 M edges) and the
"downregulates" edges swamp the informative "upregulates" edges and leave each organ node's position
governed by how many edges it has rather than by what the sample expresses. This, not the gene
weighting (Stage 1), is why liver, lung, adrenal gland, cerebellum and adipose failed on 2019.

## All 104 datasets on the 2025 'upregulates-only' build (2026-10-12)

`run_upauG_all.sh` (manifest `data/_catalog/manifest_2025upAuG.json`, three batches, 3 h 45 min), checks via
`scripts/check_build.py`, compiled by `scripts/compile_build.py --tag 2025upAuG` →
`results/compiled_2025upAuG.tsv`, `results/tissue_summary_2025upAuG.tsv`, `results/tissue_overview_2025upAuG.png`.
Each dataset is judged on its own: rank of its tissue's own gene-linked node among the build's 44 linked
Anatomy nodes (CLEAR = 1st, PARTIAL = top 7, NONE = below), plus the ontology-family check.

| | datasets |
|---|---|
| own node ranked 1st (CLEAR) | 52 |
| own node ranked 2nd–5th (PARTIAL) | 17 |
| own node below 7th (NONE) | 0 |
| tissue has no up-regulation-linked node | 35 |

- **Ranked 1st in every dataset:** liver (14 of 14), all five skeletal muscles (20 of 20, via the single
  `gastrocnemius` node), cerebrum (3), brain (2), white adipose, blood; skin 6 of 7; cerebellum 2 of 3.
- **Consistently 2nd–4th:** spleen (2nd in all 8, behind `blood`), lung (2nd–3rd in all 3, behind
  `tibial artery` / `thoracic aorta`), adrenal gland (2nd–3rd in all 4, behind `gastrocnemius`),
  hippocampus (4th in both, behind prefrontal and cingulate cortex). These are near-misses with a
  biologically related node ahead, not failures.
- **Not testable on this build** (no up-regulation edges to any node in the family): thymus, retina,
  kidney, heart, eye, bone marrow, colon, cecum, mammary gland, placenta, bone, mandible, temporal
  bone, optic nerve (35 datasets). Their ontology-family check still runs through neighbouring
  nodes: heart CLEAR in 6 of 6, colon, cecum and mandible CLEAR, kidney and retina PARTIAL, thymus NONE
  in 5 of 5 (family median 84 %).
- Analog-study controls behave like ground controls (11 CLEAR / 3 PARTIAL / 7 no node).
- Residual hub: `gastrocnemius` is in the top four of 49 datasets and `left adrenal gland cortex` of 34.
  They do not displace the tissue's own node except for adrenal gland, where gastrocnemius ranks first.

Compared with the 2019 build on the same 104 datasets (baseline-corrected, linked-node check: 6 CLEAR
out of the 31 own-node cases, connectivity Spearman −0.85), this build identifies the tissue in 69 of 69
testable datasets within the top five. The trade is coverage: 74 linked nodes in 2019 against 44 here,
with retina, kidney and heart among the tissues that lose their node.

## Cell-type route for tissues without a gene-linked Anatomy node (2026-10-12)

The 2025 graph has 2,745 CellType nodes, 59 of which carry Gene edges (33,827 `EXPRESSEDIN_GeiCT`, a
single marker-gene source, median ~400 genes each) and 846 `PARTOF_CTpA` links to Anatomy nodes.
`scripts/celltype_check.py --tag 2025upAuG` ranks all CellType nodes by baseline-corrected score in
every dataset of the upregulates-only run and reports the best-ranked cell type whose name matches the
tissue (`results/celltype_check_2025upAuG.tsv`). Ranks are out of 2,745.

| tissue (datasets) | best matching cell type | rank, median [range] |
|---|---|---|
| retina (5) | camera-type eye photoreceptor cell | 1 [1–2] |
| kidney (6) | proximal tubular cell | 1 [1–1] |
| bone marrow (2) | CD115-positive monocyte | 1 |
| placenta (1) | syncytiotrophoblast cell | 1 |
| colon (2) | smooth muscle fiber of descending colon | 1 |
| cecum (1) | Paneth cell | 2 |
| heart (6) | cardiac pacemaker cell of sinoatrial node | 2.5 [2–16] |
| mammary gland (1) | adipocyte of breast | 9 |
| thymus (5) | T cell / fetal thymocyte | 17 [3–19] |
| eye (2) | retinal rod cell | 19 [15–23] |
| optic nerve, mandible, bone, temporal bone | (chondrocyte) | 30, 47, 121, 146 |

For the tissues that do have an Anatomy node the cell-type ranking agrees with it (liver → hepatocyte
rank 1 in 14 of 14, skin → suprabasal keratinocyte, muscles → sinoatrial myocyte rank 3, cerebellum and
brain → oligodendrocyte / astrocyte), except adrenal gland, which has no steroidogenic cell type in the
59 and ranks a kidney cell at ~109.

So the cell-type route recovers retina, kidney, bone marrow, placenta, colon, cecum and heart directly
from the existing projections; thymus is weak (T cells are in the top 20 but monocytes rank above them);
eye and mammary gland are marginal; the bones and optic nerve have no cell type with gene edges and
stay untestable. Residual hub: sinoatrial-node myocytes and epicardial adipocytes appear in the top five
of about a quarter of all datasets.

## What the live graph says about the Anatomy–Gene edges (2026-10-12)

`scripts/inspect_edges.py` and `scripts/inspect_aug_coverage.py` (read-only Cypher through the SSH tunnel;
`results/anatomy_gene_edge_counts.tsv`):

| edge type | edges | Anatomy nodes | genes | properties |
|---|---|---|---|---|
| `EXPRESSES_AeG` | 1,582,757 | 92 | 19,072 | `sources` only (Bgee); blanket: blood is linked to 18,879 of 19,072 genes |
| `UPREGULATES_AuG` | 85,530 | 72 | 17,018 | `log2fc` (min 1.0, median 1.6, 90th pct 3.3), `pvalue`, `fdr`, `sources` (Bgee) |
| `DOWNREGULATES_AdG` | 493,600 | 73 | 19,014 | `log2fc`, `pvalue`, `fdr`, `sources` (Bgee) |

So a rebuild that filters the "expresses" edges by expression level is **not possible**: they carry no
value. The up-regulation edges are genuine tissue-specificity calls and they exist for 72 tissues, not
the 44 my PSEV-based detector found (it only flags nodes linked to more than half of all genes). Edge
counts per tissue range from 14 (cortex of kidney) to 6,167 (testis); heart left ventricle 484, adult
mammalian kidney 398, transverse colon 986, thoracic mammary gland 668, placenta 45, liver 113.

**Linked set re-derived from the edge list.** `check_linked.py --linked-ids` takes the set directly.
With all 72 nodes, the sparsely linked ones (heart 38 edges, brain 19, ectocervix 27) top many datasets
because a covariance over a handful of genes is noisy, so the set used is the **60 nodes with ≥ 200
up-regulation edges** (`results/aug_linked_nodes_ge200.tsv`; outputs `corrected_linked_2025upAuG60/`,
`results/compiled_2025upAuG60.tsv`, `tissue_summary_2025upAuG60.tsv`, `tissue_overview_2025upAuG60.png`).

| | 44-node set (PSEV-detected) | 60-node set (≥ 200 AuG edges) |
|---|---|---|
| own node 1st | 52 | 52 |
| own node 2nd–7th | 17 | 30 |
| own node below 7th | 0 | 2 (one colon set, brown adipose) |
| no linked node | 35 | 20 |

Newly testable with the 60-node set: heart (`heart left ventricle`, 2nd in all 6 datasets), kidney
(`adult mammalian kidney`, 1st in 5 of 6), mammary gland (5th), colon (`muscle layer of sigmoid colon`,
11th median). Still without a node: thymus, retina, eye, bone marrow, cecum, placenta (45 edges, below
the cut), bones, optic nerve; the cell-type route (previous section) covers retina, bone marrow,
placenta and cecum.

## Archive

`archive/cross_dataset/` holds the earlier analysis that compared datasets with one another
(anatomy-family heatmaps, hub-corrected shifts, specificity tests, per-dataset QQ figures, expression
ssGSEA/clustering) with its scripts, results and figures; the scripts still run from that location.
`archive/combined_check_figures/` holds the previous three-panel-in-one version of the per-dataset check.

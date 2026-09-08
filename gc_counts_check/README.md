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
run_all.sh                     projection for the six datasets + one check per dataset
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

## Archive

`archive/cross_dataset/` holds the earlier analysis that compared datasets with one another
(anatomy-family heatmaps, hub-corrected shifts, specificity tests, per-dataset QQ figures, expression
ssGSEA/clustering) with its scripts, results and figures; the scripts still run from that location.
`archive/combined_check_figures/` holds the previous three-panel-in-one version of the per-dataset check.

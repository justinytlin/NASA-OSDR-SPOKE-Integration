# NASA-OSDR-SPOKE-Integration

Tools for projecting NASA OSDR / GeneLab transcriptomes onto the
[SPOKE](https://spoke.ucsf.edu) biomedical knowledge graph, built around the
Propagated SPOKE Entry Vector (PSEV) method of
[Nelson et al., *Life* 2021](https://doi.org/10.3390/life11010042):
replicate the paper, run any new OSDR dataset through the same machinery,
control for graph artifacts, and rebuild the PSEVs from a current SPOKE.

## Which route do I need?

| I want to… | Go to | Needs |
|---|---|---|
| Run an OSDR differential-expression table through the paper's method and get ranked SPOKE nodes with a null control | [`replication/`](replication/README.md) | the 25.6 GB Zenodo PSEV zip (one-time) |
| Embed each sample of a raw counts file separately | [`replication/`](replication/README.md#counts-route-one-embedding-per-sample) | same |
| See a complete worked example, outputs included | [`examples/OSD-564/`](examples/OSD-564/README.md) | nothing |
| Ask whether a dataset's ground-control profile points at its own tissue, across many datasets | [`gc_counts_check/`](gc_counts_check/README.md) | the zip + counts files (`fetch_data.sh`) |
| Query the live SPOKE graph without any bulk data | [`spoke-api-client/`](spoke-api-client/README.md) | internet access only |
| Rebuild the gene PSEVs from a current SPOKE Neo4j instance | [`psev_build/`](psev_build/README.md) | Bolt access to a SPOKE server, ~40 GB of disk |

## Repository map

```
replication/        the PSEV pipeline: DE route, counts route, null control, anatomy check,
                    the six-study replication and its meta-analysis; figures/ for the build comparisons
examples/OSD-564/   worked example (mouse hippocampus, RR-10): rank table, null control, top findings
gc_counts_check/    ground-control tissue check across 14 OSDR datasets, one verdict per dataset
psev_build/         export a SPOKE graph over Bolt and regenerate gene PSEVs with the original algorithm
spoke-api-client/   CLI + library for the public SPOKE REST API, with neighborhood-based analyses
psev_data/          where the large downloads live (gitignored); parallel_download.sh fetches the zip
```

Every directory has its own README with the commands. Large inputs and
outputs (the PSEV zip, SPOKE node tables, built PSEVs, rank tables) are
gitignored; the READMEs say where each goes.

## Quick start

```bash
git clone https://github.com/justinytlin/NASA-OSDR-SPOKE-Integration
cd NASA-OSDR-SPOKE-Integration
git clone https://github.com/baranzini-lab/SPOKE_NASA          # node universe + paper ground truth
psev_data/parallel_download.sh                                  # 25.6 GB gene_psev.zip from Zenodo
curl -o psev_data/homologene_build68.data https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data

python3 replication/test_synthetic.py                           # smoke test, minutes, no zip needed
python3 replication/osdr_to_psev.py TABLE.csv --name OSD-NNN \
    --contrast "(Space Flight)v(Ground Control)" --filter pvalue --p 0.01
```

Requirements: Python 3.10+ with `numpy pandas scipy requests`; `matplotlib`
for figures, `scikit-learn openpyxl` for the ML and meta-analysis extras,
`neo4j` for `psev_build/`.

## The method in one paragraph

Every gene in SPOKE has a precomputed PSEV: a personalised-PageRank profile
of how strongly that gene reaches every node in the graph. A dataset's
signature genes are weighted by their log2 fold changes (or, on the counts
route, by expression), their PSEVs are summed, and every SPOKE node is
ranked, then re-ranked within its node type so that diseases, symptoms,
anatomy and pathways can each be read on their own scale. Single-study ranks
carry no significance on their own: the ground-control null control
(`replication/osdr_null_control.py`) identifies which top nodes a no-biology
split also produces, and the anatomy check uses the known source tissue as a
free positive control.

## Results at a glance

- **Replication**: the paper's headline findings reproduce on the 2019 PSEVs
  (Space Motion Sickness, vasoconstriction and cortisol-secretion processes,
  vitamin D pathways); exact numbers are unreachable because the meta-analysis
  script and two reprocessed DE tables were never published.
- **New datasets**: [`examples/OSD-564/`](examples/OSD-564/README.md) shows
  the full workflow and how to read the outputs.
- **Tissue check**: parenchymal terms of the source tissue reach the top 1–3%
  of Anatomy nodes in every dataset on the 2019 PSEVs, while broad UBERON
  families are diluted by vessel and primordium terms
  ([`gc_counts_check/README.md`](gc_counts_check/README.md)).
- **Current SPOKE**: a 2025-snapshot build (547,873 nodes) agrees with the
  2019 PSEVs at the gene, process and protein level well above chance, but
  the disease, symptom and anatomy layers were rebuilt at the source and do
  not carry over; the expression-based tissue check does not work on it and
  the reasons are documented ([`psev_build/README.md`](psev_build/README.md)).

## Data sources

| Data | Source |
|---|---|
| OSDR differential-expression and counts tables | `https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-NNN/files/` to list, `https://osdr.nasa.gov/geode-py/ws/studies/OSD-NNN/download?source=datamanager&file=...` to fetch |
| 2019 gene PSEVs (SPOKE v2, 389,297 nodes) | Zenodo [10.5281/zenodo.4404618](https://doi.org/10.5281/zenodo.4404618) |
| SPOKE v2 node table and the paper's code | [baranzini-lab/SPOKE_NASA](https://github.com/baranzini-lab/SPOKE_NASA) |
| PSEV generator (algorithm reference) | [BaranziniLab/PSEV](https://github.com/BaranziniLab/PSEV) |
| Mouse → human orthologs | NCBI HomoloGene build 68 (frozen) |
| Live SPOKE | [REST API](https://spoke.rbvi.ucsf.edu/swagger/) (public) or a Neo4j instance over Bolt (`psev_build/`) |

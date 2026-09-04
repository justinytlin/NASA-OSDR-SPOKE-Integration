#!/usr/bin/env python3
"""Route 2: raw RNA-seq counts -> per-sample full-SPOKE PSEV node ranks.

Unlike osdr_to_psev.py (which weights gene PSEVs by log2 fold change, one
embedding per contrast), this weights each gene's PSEV by that gene's
normalized expression in each individual sample, giving one embedding per
sample. NOTE: this is an extension of Nelson et al. 2021, not a replication —
the published method is fold-change-weighted.

Input: a GeneLab/OSDR unnormalized counts CSV (rows = ENSEMBL gene ids,
columns = sample names). Gene ids are mapped ENSEMBL -> mouse Entrez via a
GeneLab DE table's ENSEMBL/ENTREZID columns, then mouse -> human via frozen
HomoloGene build 68, as in the replication pipeline.

Weights per gene g, sample s:
  logcpm (default): log2(CPM[g,s] + 1)
  zscore:           log2(CPM+1) z-scored per gene across the kept samples
                    (the spoke_embeddings_ml.py convention, minus the
                    up/down split; signed weights are fine downstream)

Output matches the replication rank-table shape (node_info columns + one
overall-rank and one Rank_by_type_ column per sample), so welch_nodes.py /
welch_meta.py-style tooling works on it.

Usage:
  counts_to_psev.py COUNTS.csv --name OSD-564-GC --samples _GC_
      [--gene-map DE_TABLE.csv]   ENSEMBL->ENTREZID source (default: OSD-564 DE table)
      [--weight logcpm|zscore]    per-gene weights (default logcpm)
      [--min-cpm X]               drop genes with CPM < X in every kept sample (default 1.0)
      [--out DIR] [--zip PSEV_ZIP] [--spoke SPOKE_DIR] [--homologene FILE]
      [-b 0.1] [--prep-only] [--limit-groups N (debug)]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(HERE))
from psev_pipeline import compute_rank_table, load_homologene_mouse_to_human


def load_counts(path, sample_filters):
    df = pd.read_csv(path, header=0, index_col=0)
    if sample_filters:
        keep = [c for c in df.columns if any(s in c for s in sample_filters)]
        if not keep:
            sys.exit(f"no sample columns matched {sample_filters}; "
                     f"available: {list(df.columns)}")
        df = df[keep]
    return df


def load_ensembl_to_entrez(gene_map_path):
    with open(gene_map_path) as fh:
        first = fh.readline()
    sep = "\t" if first.count("\t") > first.count(",") else ","
    m = pd.read_csv(gene_map_path, sep=sep, usecols=["ENSEMBL", "ENTREZID"],
                    low_memory=False)
    return m[m.ENTREZID.notna()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts", help="unnormalized counts CSV (ENSEMBL x samples)")
    ap.add_argument("--name", default=None,
                    help="label for the output file (default: from filename)")
    ap.add_argument("--samples", action="append", default=None,
                    help="substring filter on sample columns; repeatable")
    ap.add_argument("--gene-map",
                    default=str(BASE / "psev_data/OSD-564_rna_seq_differential_expression.csv"),
                    help="table with ENSEMBL and ENTREZID columns")
    ap.add_argument("--weight", choices=["logcpm", "zscore"], default="logcpm")
    ap.add_argument("--min-cpm", type=float, default=1.0,
                    help="keep genes with CPM >= this in at least one kept sample")
    ap.add_argument("--out", default=str(HERE / "psev_out"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    ap.add_argument("--prep-only", action="store_true")
    ap.add_argument("--limit-groups", type=int, default=None)
    args = ap.parse_args()

    name = args.name or Path(args.counts).stem.split("_")[0]
    save_str = "_".join(str(args.b).split("."))
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)

    node_info_df = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    spoke_genes = set(node_info_df[node_info_df.Node_Type == "Gene"].Node.values)

    counts = load_counts(args.counts, args.samples)
    sample_cols = list(counts.columns)
    print(f"[{name}] {Path(args.counts).name}: {len(counts)} genes, "
          f"{len(sample_cols)} samples kept: {sample_cols}")

    # CPM -> log2(CPM+1); drop genes never reaching --min-cpm
    cpm = counts / counts.sum(axis=0) * 1e6
    expressed = (cpm >= args.min_cpm).any(axis=1)
    weights = np.log2(cpm[expressed] + 1.0)
    print(f"  genes with CPM>={args.min_cpm} in >=1 sample: {len(weights)}")

    # ENSEMBL -> mouse Entrez -> human Entrez (SPOKE Node)
    e2e = load_ensembl_to_entrez(args.gene_map)
    e2e.loc[:, "ENTREZID"] = e2e.ENTREZID.astype(str).str.split("|")
    e2e = e2e.explode("ENTREZID")
    e2e = e2e[e2e.ENTREZID.str.fullmatch(r"\d+(\.0)?")]
    e2e.loc[:, "ENTREZID"] = e2e.ENTREZID.astype(float).astype(int)
    m2h = load_homologene_mouse_to_human(args.homologene)
    e2e = e2e.merge(m2h[["ENTREZID", "Human_EntrezGene ID"]], on="ENTREZID")
    e2e.loc[:, "Node"] = e2e["Human_EntrezGene ID"].astype(int).astype(str)
    e2e = e2e[e2e.Node.isin(spoke_genes)][["ENSEMBL", "Node"]].drop_duplicates()

    kept = weights.reset_index().rename(columns={weights.index.name or "index": "ENSEMBL"})
    kept = kept.merge(e2e, on="ENSEMBL")
    kept = kept.groupby("Node")[sample_cols].mean().reset_index()
    print(f"  mapped SPOKE genes: {len(kept)}")

    if args.weight == "zscore":
        vals = kept[sample_cols].values
        mu = vals.mean(axis=1, keepdims=True)
        sd = vals.std(axis=1, keepdims=True)
        sd[sd == 0] = np.nan
        kept.loc[:, sample_cols] = np.nan_to_num((vals - mu) / sd)
        kept = kept[np.abs(kept[sample_cols].values).max(axis=1) > 0]
        print(f"  z-scored across samples; nonconstant genes: {len(kept)}")

    if not len(kept):
        sys.exit("no genes left")
    if args.prep_only:
        return

    result_df = compute_rank_table(kept, sample_cols, node_info_df, args.zip,
                                   b=args.b, limit_groups=args.limit_groups)
    out = out_dir / f"{name}_counts_ranks_and_rank_by_type_for_meta{save_str}.tsv"
    result_df.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

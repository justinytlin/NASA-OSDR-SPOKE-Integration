#!/usr/bin/env python3
"""Bridge: OSDR/GeneLab differential-expression table -> full-SPOKE PSEV node ranks.

Feeds the DE tables the spoke-api-client examples already consume (CSV from
https://visualization.osdr.nasa.gov/biodata/api/v2/... or the geode-py study
download, or the repo's GeneLab TSVs) through the Nelson-2021 replication
machinery: mouse->human orthologs from frozen HomoloGene build 68, |log2FC|-
weighted sum of precomputed SPOKE v2 gene PSEVs, overall + per-node-type ranks
over all 389,297 nodes.

Output matches the replication's <name>_ranks_and_rank_by_type_for_meta<b>.tsv
shape, so welch_meta.py-style meta-analysis works downstream unchanged.

Usage:
  osdr_to_psev.py TABLE.csv --name OSD-352
      [--human]                   input genes are human; use ENTREZID directly
      [--contrast SUBSTR ...]     keep only Log2fc_ columns containing SUBSTR
      [--filter none|pvalue|same-direction]   gene selection (default none)
      [--p 0.01 | --fdr 0.05]     cutoffs for --filter pvalue
      [--out DIR] [--zip PSEV_ZIP] [--spoke SPOKE_DIR] [--homologene FILE]
      [-b 0.1] [--prep-only] [--limit-groups N (debug)]
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(HERE))
from psev_pipeline import (compute_rank_table, filter_same_direction,
                           load_homologene_mouse_to_human)


def load_de_table(path):
    """OSDR CSV or GeneLab TSV; strips 'a/b/COL' header prefixes like the
    spoke-api-client loaders do."""
    with open(path) as fh:
        first = fh.readline()
    sep = "\t" if first.count("\t") > first.count(",") else ","
    df = pd.read_csv(path, sep=sep, header=0, index_col=False, low_memory=False)
    df.columns = [c.split("/")[-1] for c in df.columns]
    return df


def map_to_spoke_genes(exp_df, spoke_genes, homologene_path, human=False):
    """Attach a SPOKE gene-node id ('Node', human Entrez as str) to each row."""
    exp_df = exp_df[exp_df.ENTREZID.notna()].copy()
    # newer GLbulkRNAseq tables carry multi-mapped ids like "54204|100043580"
    exp_df.loc[:, "ENTREZID"] = exp_df.ENTREZID.astype(str).str.split("|")
    exp_df = exp_df.explode("ENTREZID")
    exp_df = exp_df[exp_df.ENTREZID.str.fullmatch(r"\d+(\.0)?")]
    exp_df.loc[:, "ENTREZID"] = exp_df.ENTREZID.astype(float).astype(int)
    if human:
        exp_df.loc[:, "Node"] = exp_df.ENTREZID.astype(str)
    else:
        m2h = load_homologene_mouse_to_human(homologene_path)
        exp_df = exp_df.merge(
            m2h[["ENTREZID", "Human_EntrezGene ID"]], on="ENTREZID")
        exp_df.loc[:, "Node"] = (
            exp_df["Human_EntrezGene ID"].astype(int).astype(str))
    return exp_df[exp_df.Node.isin(spoke_genes)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", help="OSDR/GeneLab differential expression table")
    ap.add_argument("--name", default=None,
                    help="study label for the output file (default: from filename)")
    ap.add_argument("--human", action="store_true",
                    help="genes are already human (e.g. OSD-871); skip orthologs")
    ap.add_argument("--contrast", action="append", default=None,
                    help="substring filter on Log2fc_ columns; repeatable")
    ap.add_argument("--filter", choices=["none", "pvalue", "same-direction"],
                    default="none")
    ap.add_argument("--p", type=float, default=0.01)
    ap.add_argument("--fdr", type=float, default=None)
    ap.add_argument("--out", default=str(HERE / "psev_out"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    ap.add_argument("--prep-only", action="store_true")
    ap.add_argument("--limit-groups", type=int, default=None)
    args = ap.parse_args()

    name = args.name or Path(args.table).stem.split("_")[0]
    save_str = "_".join(str(args.b).split("."))
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)

    node_info_df = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    spoke_genes = set(node_info_df[node_info_df.Node_Type == "Gene"].Node.values)

    exp_df = load_de_table(args.table)
    fc_cols = [c for c in exp_df.columns if c.startswith("Log2fc_")]
    if args.contrast:
        fc_cols = [c for c in fc_cols if any(s in c for s in args.contrast)]
    if not fc_cols:
        sys.exit("no Log2fc_ columns matched")
    stat_cols = [c for c in exp_df.columns
                 if re.match(r"P\.value_|Adj\.p\.value_|LRT\.p\.value", c)]
    print(f"[{name}] {Path(args.table).name}: {len(exp_df)} rows, "
          f"{len(fc_cols)} contrasts selected")

    exp_df = map_to_spoke_genes(exp_df[["ENTREZID"] + fc_cols + stat_cols],
                                spoke_genes, args.homologene, human=args.human)
    # mean for genes seen more than once; drop rows with no usable fc
    exp_df = exp_df.groupby("Node").mean(numeric_only=True).reset_index()
    exp_df[fc_cols] = exp_df[fc_cols].fillna(0.0)
    exp_df = exp_df[np.abs(exp_df[fc_cols].values).max(axis=1) > 0]
    print(f"  mapped SPOKE genes: {len(exp_df)}")

    if args.filter == "pvalue":
        prefix = "Adj.p.value_" if args.fdr is not None else "P.value_"
        cutoff = args.fdr if args.fdr is not None else args.p
        pcols = [prefix + c[len("Log2fc_"):] for c in fc_cols]
        pcols = [c for c in pcols if c in exp_df.columns]
        if not pcols:
            sys.exit(f"no {prefix} columns for the selected contrasts")
        keep = np.nanmin(exp_df[pcols].values, axis=1) < cutoff
        kept = exp_df[keep]
        print(f"  {prefix}<{cutoff} in any contrast: {len(kept)} genes kept")
    elif args.filter == "same-direction":
        exp_df.loc[:, "max_fc"] = np.max(exp_df[fc_cols].values, axis=1)
        exp_df = filter_same_direction(exp_df, np.array(fc_cols))
        kept = exp_df[exp_df.Space_over_Ground_Basal_same_sign == True]
        print(f"  same-direction genes kept: {len(kept)}")
    else:
        kept = exp_df
        print(f"  no gene filter: {len(kept)} genes kept")
    if not len(kept):
        sys.exit("no genes left after filtering")
    if args.prep_only:
        return

    result_df = compute_rank_table(kept, fc_cols, node_info_df, args.zip,
                                   b=args.b, limit_groups=args.limit_groups)
    out = out_dir / f"{name}_ranks_and_rank_by_type_for_meta{save_str}.tsv"
    result_df.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

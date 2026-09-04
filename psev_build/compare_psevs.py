#!/usr/bin/env python3
"""Compare gene PSEVs from two stores (zip or directory) over their shared nodes.

For each gene present in both stores: Spearman correlation of the two PSEVs
restricted to nodes shared by both node universes (matched on
(Node identifier, Node_Type)), plus overlap of the top-K shared nodes.

Usage:
  compare_psevs.py --a STORE_A --a-nodes node_info_A.tsv
                   --b STORE_B --b-nodes node_info_B.tsv
                   [--genes 1,10,100 | --n 20] [--top 50] [--jump 0.1]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "replication"))
from psev_pipeline import PsevStore, load_gene_indices_df  # noqa: E402


def load(store_path, nodes_path):
    nodes = pd.read_csv(nodes_path, sep="\t", dtype={"Node": str})
    store = PsevStore(store_path)
    gi = load_gene_indices_df(store, nodes.Node.values)
    gi = gi.set_index("Node")
    return nodes, store, gi


def get_row(store, gi, gene, save_str, cache):
    r = int(gi.loc[gene, "Round"]); ri = int(gi.loc[gene, "round_index"])
    if r not in cache:
        cache.clear()
        cache[r] = store.load_group(r, save_str)
    return np.asarray(cache[r][ri], dtype=np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True); ap.add_argument("--a-nodes", required=True)
    ap.add_argument("--b", required=True); ap.add_argument("--b-nodes", required=True)
    ap.add_argument("--genes", default=None); ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--top", type=int, default=50); ap.add_argument("--jump", type=float, default=0.1)
    args = ap.parse_args()
    save_str = "_".join(str(args.jump).split("."))

    na, sa, ga = load(args.a, args.a_nodes)
    nb, sb, gb = load(args.b, args.b_nodes)
    key = lambda df: df.Node.astype(str) + "|" + df.Node_Type.astype(str)
    ka, kb = key(na), key(nb)
    shared = np.intersect1d(ka.values, kb.values)
    ia = pd.Series(np.arange(len(na)), index=ka.values).loc[shared].values
    ib = pd.Series(np.arange(len(nb)), index=kb.values).loc[shared].values
    print(f"universe A {len(na):,}  B {len(nb):,}  shared {len(shared):,}")
    types = na.Node_Type.values[ia]

    genes = args.genes.split(",") if args.genes else \
        [g for g in ga.index if g in gb.index][:args.n]
    ca, cb = {}, {}
    rows = []
    for g in genes:
        if g not in ga.index or g not in gb.index:
            print(f"{g}: missing in one store"); continue
        va = get_row(sa, ga, g, save_str, ca)[ia]
        vb = get_row(sb, gb, g, save_str, cb)[ib]
        rho = spearmanr(va, vb).correlation
        ta = set(np.argsort(-va)[:args.top]); tb = set(np.argsort(-vb)[:args.top])
        rows.append((g, rho, len(ta & tb) / args.top))
        top_a = [f"{na.Node_Name.values[ia][i]}({types[i]})" for i in np.argsort(-va)[:5]]
        top_b = [f"{nb.Node_Name.values[ib][i]}({types[i]})" for i in np.argsort(-vb)[:5]]
        print(f"{g}: rho={rho:.3f} top{args.top} overlap={rows[-1][2]:.2f}\n   A: {top_a}\n   B: {top_b}")
    df = pd.DataFrame(rows, columns=["gene", "spearman", f"top{args.top}_overlap"])
    print(df.describe().loc[["mean", "50%", "min", "max"]].to_string())


if __name__ == "__main__":
    main()

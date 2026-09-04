#!/usr/bin/env python3
"""Build gene PSEVs (Propagated SPOKE Entry Vectors) from an exported SPOKE graph.

Reimplements BaranziniLab/PSEV make_psevs_by_node_type.py (the generator behind
the Zenodo gene_psev.zip used by Nelson et al. 2021) on a current SPOKE export:

  * adjacency A: binary, symmetric (undirected), parallel edges collapsed,
    self loops dropped -- exactly what the original got from neo4j_edges.tsv
  * transition matrix: A row-normalised by degree (dangling nodes -> zero row)
  * restart: one-hot on the gene (jump probability b = 0.1)
  * power iteration from the uniform vector, stop when the L1 change <= a
    (0.001) or after 40 iterations, then renormalise to sum 1
    (--converge instead iterates to L1 change <= 1e-12 / 1000 iterations)

Output layout mirrors gene_psev.zip so replication/psev_pipeline.py can read it
directly (pass the output directory as --zip):

  <out>/gene_psevs/gene_group_<g>.tsv                 one column node_2_index
  <out>/gene_psevs/raw_psev_0_1_gene_group_<g>_sparse.npy   (genes_in_group, n_nodes)

Usage:
  make_gene_psevs.py --spoke DIR [--out DIR] [--group-size 1000] [--dtype float32]
      [--node-types A,B,...]        restrict the node universe (edges to other
                                    types are dropped); writes node_info.tsv for
                                    the restricted universe into --out
      [--exclude-rel-types X,Y]     ignore some relationship types
      [--prune-types Compound]      keep those types only where they touch another type
      [--drop-isolated]             drop nodes without edges in the universe
      [--genes FILE | --gene-limit N]   subset of Gene identifiers (one per line)
      [--workers 8] [--block 16] [-b 0.1] [-a 0.001] [--max-iter 40] [--converge]
      [--dry-run]                   print sizes and exit
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent

# shared (fork-inherited) state for workers
_P = None          # csr transition matrix (row-stochastic on the symmetric graph)
_CFG = {}


def build_graph(spoke_dir, node_types=None, exclude_rel=None, prune_types=None, drop_isolated=False):
    nodes = pd.read_csv(spoke_dir / "node_info.tsv", sep="\t", dtype={"Node": str})
    iids = np.load(spoke_dir / "node_internal_ids.npy")
    assert len(iids) == len(nodes)
    if node_types:
        keep = nodes.Node_Type.isin(node_types).values
        nodes = nodes[keep].reset_index(drop=True)
        iids = iids[keep]
        nodes["Node_Index"] = np.arange(len(nodes))
    n = len(nodes)
    # internal id -> index lookup (internal ids are dense-ish ints)
    lut = np.full(int(iids.max()) + 1, -1, dtype=np.int64)
    lut[iids] = np.arange(n)

    edir = spoke_dir / "edges"
    files = sorted(edir.glob("*.npy"))
    if exclude_rel:
        files = [f for f in files if f.stem not in exclude_rel]
    src, dst = [], []
    per_type = {}
    for f in files:
        e = np.load(f)
        if len(e) == 0:
            continue
        ok = (e < len(lut)).all(axis=1)
        e = e[ok]
        s, t = lut[e[:, 0]], lut[e[:, 1]]
        m = (s >= 0) & (t >= 0) & (s != t)
        per_type[f.stem] = int(m.sum())
        src.append(s[m]); dst.append(t[m])
    src = np.concatenate(src); dst = np.concatenate(dst)

    def adjacency(src, dst, n):
        # symmetric binary adjacency; duplicates collapse via csr + data=1
        rows = np.concatenate([src, dst]); cols = np.concatenate([dst, src])
        A = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, n))
        A.data[:] = 1.0
        A.sum_duplicates()
        return A

    A = adjacency(src, dst, n)
    pruned = {}
    keep = np.ones(n, dtype=bool)
    if prune_types:
        # drop nodes of these types that have no edge to a node of another type
        # (e.g. Compound nodes hanging only off the ChEBI is-a tree or isolated)
        types = nodes.Node_Type.values
        for t in prune_types:
            is_t = types == t
            other = (~is_t).astype(np.float32)
            n_other = np.asarray(A @ other).ravel()   # neighbours of another type
            drop = is_t & (n_other == 0)
            pruned[t] = int(drop.sum())
            keep &= ~drop
    if drop_isolated:
        deg0 = np.asarray(A.sum(axis=1)).ravel()
        iso = deg0 == 0
        pruned["isolated_any_type"] = int((iso & keep).sum())
        keep &= ~iso
    if not keep.all():
        remap = np.full(n, -1, dtype=np.int64); remap[keep] = np.arange(keep.sum())
        m = keep[src] & keep[dst]
        src, dst = remap[src[m]], remap[dst[m]]
        nodes = nodes[keep].reset_index(drop=True); nodes["Node_Index"] = np.arange(len(nodes))
        n = len(nodes)
        A = adjacency(src, dst, n)
        print(f"pruned: {pruned}")
    deg = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(deg); nz = deg > 0; inv[nz] = 1.0 / deg[nz]
    P = sp.diags(inv) @ A.astype(np.float64)          # row-stochastic, float64
    P = P.tocsr(); P.sort_indices()
    info = {"n_nodes": int(n), "undirected_edges": int(A.nnz // 2),
            "isolated_nodes": int((~nz).sum()), "edges_used_per_type": per_type,
            "pruned": pruned,
            "nodes_per_type": nodes.Node_Type.value_counts().to_dict()}
    return nodes, P, deg, info


def ppr_block(restart_idx):
    """Personalised PageRank for a block of restart nodes (columns)."""
    P, cfg = _P, _CFG
    n = P.shape[0]; k = len(restart_idx)
    b, a, max_iter, dtype = cfg["b"], cfg["a"], cfg["max_iter"], cfg["dtype"]
    R = np.zeros((n, k), dtype=dtype)
    R[restart_idx, np.arange(k)] = b
    X = np.full((n, k), 1.0 / n, dtype=dtype)
    PT = cfg["PT"]
    active = np.ones(k, dtype=bool)
    it = 0
    while active.any() and it < max_iter:
        # original: new = ((b*r + (1-b)*M).T) @ v  ->  b*r*sum(v) + (1-b)*M.T@v
        # sum(v) drifts below 1 when mass leaks through dangling nodes
        Xn = R * X.sum(axis=0, keepdims=True) + (1.0 - b) * (PT @ X)
        diff = np.abs(Xn - X).sum(axis=0)
        X = Xn
        active = diff > a
        it += 1
    X /= X.sum(axis=0, keepdims=True)
    return X.T.astype(cfg["out_dtype"])


def _init(PT, cfg):
    global _P, _CFG
    _P = PT  # only shape used
    _CFG = dict(cfg, PT=PT)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spoke", required=True, help="export dir from export_spoke.py")
    ap.add_argument("--out", default=None, help="default: <spoke>/gene_psev")
    ap.add_argument("--node-types", default=None)
    ap.add_argument("--exclude-rel-types", default=None)
    ap.add_argument("--prune-types", default=None,
                    help="comma-separated node types to keep only when they touch another type (e.g. Compound)")
    ap.add_argument("--drop-isolated", action="store_true", help="drop nodes with no edges in the universe")
    ap.add_argument("--genes", default=None, help="file of Gene identifiers to compute")
    ap.add_argument("--gene-limit", type=int, default=None)
    ap.add_argument("--group-size", type=int, default=1000)
    ap.add_argument("--dtype", default="float32", choices=["float32", "float64"],
                    help="storage dtype (original zip is float64; float32 halves disk; computation is always float64)")
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1))
    ap.add_argument("--block", type=int, default=16, help="restart nodes per worker task")
    ap.add_argument("-b", type=float, default=0.1, help="random-jump probability")
    ap.add_argument("-a", type=float, default=0.001, help="L1 stopping threshold")
    ap.add_argument("--max-iter", type=int, default=40)
    ap.add_argument("--converge", action="store_true", help="a=1e-12, max-iter=1000")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    spoke = Path(args.spoke)
    out = Path(args.out) if args.out else spoke / "gene_psev"
    split = lambda s: set(x.strip() for x in s.split(",") if x.strip()) if s else None
    if args.converge:
        args.a, args.max_iter = 1e-12, 1000

    t0 = time.time()
    nodes, P, deg, info = build_graph(spoke, split(args.node_types), split(args.exclude_rel_types),
                                      split(args.prune_types), args.drop_isolated)
    n = len(nodes)
    print(f"graph: {n:,} nodes, {info['undirected_edges']:,} undirected edges, "
          f"{info['isolated_nodes']:,} isolated ({time.time() - t0:.0f}s)")
    print("nodes per type:", info["nodes_per_type"])
    print("edges used per type:", info["edges_used_per_type"])

    genes = nodes[nodes.Node_Type == "Gene"]
    if args.genes:
        want = set(Path(args.genes).read_text().split())
        genes = genes[genes.Node.isin(want)]
    if args.gene_limit:
        genes = genes.iloc[:args.gene_limit]
    gene_idx = genes.Node_Index.values.astype(np.int64)
    n_groups = int(np.ceil(len(gene_idx) / args.group_size))
    itemsize = np.dtype(args.dtype).itemsize
    total_gb = len(gene_idx) * n * itemsize / 1e9
    print(f"genes: {len(gene_idx):,} -> {n_groups} groups of {args.group_size}; "
          f"output ~{total_gb:,.1f} GB as {args.dtype}")
    print(f"restart b={args.b}, stop a={args.a}, max_iter={args.max_iter}, "
          f"workers={args.workers}, block={args.block}")
    if args.dry_run:
        return

    out.mkdir(parents=True, exist_ok=True)
    gdir = out / "gene_psevs"; gdir.mkdir(exist_ok=True)
    nodes.to_csv(out / "node_info.tsv", sep="\t", index=False)
    info.update({"b": args.b, "a": args.a, "max_iter": args.max_iter, "dtype": args.dtype,
                 "n_genes": int(len(gene_idx)), "group_size": args.group_size,
                 "node_types": sorted(nodes.Node_Type.unique().tolist()),
                 "source": str(spoke)})
    (out / "psev_build_info.json").write_text(json.dumps(info, indent=2))

    PT = P.T.tocsr()
    del P
    # iterate and normalise in float64 regardless of the storage dtype: a float32
    # column sum over ~5e5 nodes is only accurate to ~1e-3, which showed up as
    # row sums of 1.0017 in the first build
    cfg = {"b": args.b, "a": args.a, "max_iter": args.max_iter,
           "dtype": np.float64, "out_dtype": np.dtype(args.dtype)}
    save_str = "_".join(str(args.b).split("."))

    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init, initargs=(PT, cfg)) as pool:
        for g in range(n_groups):
            tsv = gdir / f"gene_group_{g}.tsv"
            npy = gdir / f"raw_psev_{save_str}_gene_group_{g}_sparse.npy"
            if npy.exists() and tsv.exists():
                print(f"group {g}: exists, skipping")
                continue
            gi = gene_idx[g * args.group_size:(g + 1) * args.group_size]
            tg = time.time()
            blocks = [gi[i:i + args.block] for i in range(0, len(gi), args.block)]
            mat = np.empty((len(gi), n), dtype=args.dtype)
            pos = 0
            for res in pool.imap(ppr_block, blocks):
                mat[pos:pos + len(res)] = res
                pos += len(res)
            tmp = npy.with_suffix(".tmp.npy")
            np.save(tmp, mat)
            tmp.rename(npy)
            pd.DataFrame({"node_2_index": gi}).to_csv(tsv, sep="\t", index=False)
            del mat
            print(f"group {g}/{n_groups - 1}: {len(gi)} genes in {time.time() - tg:.0f}s "
                  f"(elapsed {(time.time() - t0) / 60:.1f} min)", flush=True)
    print(f"done: {out}")


if __name__ == "__main__":
    main()

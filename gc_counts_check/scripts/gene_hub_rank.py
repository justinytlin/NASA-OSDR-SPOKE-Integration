#!/usr/bin/env python3
"""Per-gene hub-rank covariate for degree-matched permutation nulls.

For every gene PSEV (all 19,567 genes, z-scored per node over all genes, ranked within
the gene exactly as in the projection) record
  hub_rank  = mean within-gene rank (0-1 scale) of the 74 gene-linked Anatomy nodes
  n_high_z  = number of nodes with z > 3  (how many strong associations the gene has)
A gene that assigns high ranks to the hub nodes pushes them up whenever it carries
weight; permuting weights among genes with similar hub_rank removes that effect.
Output psev_out/_gene_hub_rank.tsv.  ~10 min.
"""
import sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent.parent; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "replication")); sys.path.insert(0, str(ROOT / "scripts"))
from psev_pipeline import PsevStore, load_gene_indices_df
from gc_counts_to_psev import rank_rows

ni = pd.read_csv(BASE / "SPOKE_NASA/spoke_v_2/node_info.tsv", sep="\t"); n_nodes = len(ni)
A = np.load(ROOT / "psev_out/_anatomy_psev.npz")
link_idx = A["anatomy_node_index"][(A["values"] > 1e-4).mean(axis=0) > 0.5]
t0 = time.time(); s = np.zeros(n_nodes); ss = np.zeros(n_nodes); out = {}
with PsevStore(BASE / "psev_data/gene_psev.zip") as zf, ThreadPoolExecutor(7) as pool:
    gi = load_gene_indices_df(zf, ni.Node.values); rounds = sorted(gi.Round.unique()); total = len(gi)
    for r in rounds:
        mat = zf.load_group(r, "0_1")
        for i in range(0, mat.shape[0], 128):
            sub = mat[i:i + 128].astype(np.float64); s += sub.sum(axis=0); ss += (sub ** 2).sum(axis=0)
        del mat; print(f"pass1 {r} ({time.time() - t0:.0f}s)", flush=True)
    avg = s / total; std = np.sqrt(ss / total - avg ** 2); std[std == 0] = np.nan

    def work(mat, idx):
        z = np.nan_to_num((mat[idx].astype(np.float64) - avg) / std)
        rk = rank_rows(z) / (n_nodes - 1)
        return idx, rk[:, link_idx].mean(axis=1), (z > 3).sum(axis=1)

    for r in rounds:
        mat = zf.load_group(r, "0_1"); g = gi[gi.Round == r].sort_values("round_index")
        futs = [pool.submit(work, mat, np.arange(i, min(i + 64, mat.shape[0]))) for i in range(0, mat.shape[0], 64)]
        for f in futs:
            idx, h, nz = f.result()
            for j, hv, nv in zip(idx, h, nz):
                out[(r, int(j))] = (hv, int(nv))
        del mat; print(f"pass2 {r} ({time.time() - t0:.0f}s)", flush=True)
gi["hub_rank"] = [out[(r, j)][0] for r, j in zip(gi.Round, gi.round_index)]
gi["n_high_z"] = [out[(r, j)][1] for r, j in zip(gi.Round, gi.round_index)]
gi[["Node", "hub_rank", "n_high_z"]].to_csv(ROOT / "psev_out/_gene_hub_rank.tsv", sep="\t", index=False)
print("wrote psev_out/_gene_hub_rank.tsv", len(gi))

#!/usr/bin/env python3
"""Integrity check of a gene-PSEV build directory: every group file loads, its
shape matches its gene_group TSV, sampled rows are finite, non-negative, sum to
1 (within float32 rounding) and put >= b on the restart gene.

Usage: check_build.py BUILD_DIR [--rows-per-group 10] [--jump 0.1]
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "replication"))
from psev_pipeline import PsevStore, load_gene_indices_df  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("build"); ap.add_argument("--rows-per-group", type=int, default=10)
ap.add_argument("--jump", type=float, default=0.1)
a = ap.parse_args()
B = Path(a.build); save_str = "_".join(str(a.jump).split("."))
nodes = pd.read_csv(B / "node_info.tsv", sep="\t", dtype={"Node": str}); n = len(nodes)
st = PsevStore(B); gi = load_gene_indices_df(st, nodes.Node.values)
print(f"universe {n:,}  genes with PSEVs {len(gi):,}  groups {len(st.groups())}")
bad = 0; t0 = time.time()
for r in st.groups():
    m = st.load_group(r, save_str); g = gi[gi.Round == r]
    if m.shape != (len(g), n):
        print(f"group {r}: SHAPE MISMATCH {m.shape} vs ({len(g)}, {n})"); bad += 1; continue
    rows = np.unique(np.linspace(0, len(g) - 1, a.rows_per_group).astype(int))
    sub = np.asarray(m[rows], dtype=np.float64)
    sums = sub.sum(axis=1); selfv = sub[np.arange(len(rows)), g.node_2_index.values[rows]]
    tol = 1e-5 if m.dtype == np.float32 else 1e-9
    ok = (np.isfinite(sub).all() and np.allclose(sums, 1, atol=tol)
          and (selfv >= a.jump * 0.999).all() and (sub.min() >= 0))
    bad += (not ok)
    print(f"group {r:2d}: {m.shape} {m.dtype} sums {sums.min():.7f}-{sums.max():.7f} "
          f"self {selfv.min():.4f}-{selfv.max():.4f} min {sub.min():.1e} {'OK' if ok else 'BAD'}")
print(f"groups with problems: {bad} ({time.time() - t0:.0f}s)")
sys.exit(1 if bad else 0)

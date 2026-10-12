#!/usr/bin/env python3
"""Extract the Anatomy columns of every gene PSEV (2019 SPOKE) into one matrix.

Output psev_out/_anatomy_psev.npz: values float32 [n_genes x 13,257 Anatomy nodes],
gene_node_index (row -> SPOKE node index of the gene), self_value (the gene's PSEV at itself).
Used to find which Anatomy nodes carry gene-expression links (the edge list of the 2019
graph is not published with the PSEVs) and for no-propagation controls. ~5 min, 1 GB.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent.parent; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "replication"))
from psev_pipeline import PsevStore, load_gene_indices_df

ap = argparse.ArgumentParser()
ap.add_argument("--build", default=str(BASE / "psev_data/gene_psev.zip"), help="PSEV zip or build directory")
ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"), help="directory with node_info.tsv")
ap.add_argument("--out", default=str(ROOT / "psev_out/_anatomy_psev.npz"))
args = ap.parse_args()
ni = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
ana = np.where(ni.Node_Type.values == "Anatomy")[0]
vals, gidx, selfv = [], [], []
t0 = time.time()
with PsevStore(args.build) as zf:
    gi = load_gene_indices_df(zf, ni.Node.values)
    for r in sorted(gi.Round.unique()):
        g = gi[gi.Round == r].sort_values("round_index")
        mat = zf.load_group(r, "0_1")
        vals.append(mat[:, ana].astype(np.float32)); gidx.append(g.node_2_index.values)
        selfv.append(mat[np.arange(len(g)), g.node_2_index.values].astype(np.float32))
        del mat; print(f"group {r} done ({time.time() - t0:.0f}s)", flush=True)
np.savez(args.out, values=np.vstack(vals), gene_node_index=np.concatenate(gidx),
         self_value=np.concatenate(selfv), anatomy_node_index=ana)
print("wrote", args.out, np.vstack(vals).shape)

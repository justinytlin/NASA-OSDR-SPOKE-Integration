#!/usr/bin/env python3
"""Per-sample SPOKE embeddings from OSD-871 + leave-one-out classification.

Method (PSEV-style, after Nelson et al. 2019):
  1. Subgraph: union of 1-hop SPOKE neighborhoods of the signature genes
     (fetched by fetch_subgraph.py).
  2. Per sample: z-score each signature gene's expression across samples;
     split into up (z>0) and down (z<0) entry weights on the gene nodes.
  3. Personalized PageRank (restart prob. alpha) with each half as the
     restart distribution -> two vectors over all subgraph nodes; concatenate.
  4. Leave-one-out logistic regression on the embeddings for:
     flight vs ground, PD vs healthy donor, with/without microglia.
  5. PCA of the embedding space -> pca_coords.csv (+ scatter PNG).

Usage:
  spoke_embeddings_ml.py --subgraph osd871_ml/subgraph.jsonl \
      --de TABLE.csv --samples SAMPLES.csv --out osd871_ml
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import sparse

ALPHA = 0.15   # restart probability
TOL, MAX_ITER = 1e-10, 200


def personalized_pagerank(P, restart, alpha=ALPHA):
    """Power iteration for PPR. P: column-stochastic sparse matrix."""
    r = restart / restart.sum()
    x = r.copy()
    for _ in range(MAX_ITER):
        x_new = (1 - alpha) * (P @ x) + alpha * r
        if np.abs(x_new - x).sum() < TOL:
            return x_new
        x = x_new
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subgraph", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(exist_ok=True)

    # ── build the graph ────────────────────────────────────────────────
    node_index, node_meta, edge_set = {}, [], set()
    gene_node = {}
    for line in open(args.subgraph):
        rec = json.loads(line)
        for n in rec["nodes"]:
            if n["id"] not in node_index:
                node_index[n["id"]] = len(node_meta)
                node_meta.append(n)
            if n["type"] == "Gene" and n["name"] == rec["gene"]:
                gene_node[rec["gene"]] = node_index[n["id"]]
        for e in rec["edges"]:
            edge_set.add((e["s"], e["t"]))
    N = len(node_meta)
    rows, cols = [], []
    for s, t in edge_set:
        si, ti = node_index.get(s), node_index.get(t)
        if si is None or ti is None:
            continue
        rows += [si, ti]        # undirected
        cols += [ti, si]
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N))
    deg = np.asarray(A.sum(axis=0)).ravel()
    deg[deg == 0] = 1
    P = A @ sparse.diags(1.0 / deg)   # column-stochastic
    print(f"subgraph: {N} nodes, {len(edge_set)} edges, "
          f"{len(gene_node)} signature genes anchored")

    # ── per-sample expression z-scores for signature genes ─────────────
    r = csv.reader(open(args.de))
    hdr = next(r)
    ix = {h: i for i, h in enumerate(hdr)}
    sample_cols = [h for h in hdr if h.startswith("GSM")]
    expr = {}
    for row in r:
        sym = row[ix["SYMBOL"]]
        if sym in gene_node and sym not in expr:
            try:
                expr[sym] = np.array([float(row[ix[c]]) for c in sample_cols])
            except ValueError:
                continue
    genes = sorted(expr)
    X = np.log2(np.vstack([expr[g] for g in genes]) + 1.0)
    Z = (X - X.mean(axis=1, keepdims=True)) / (X.std(axis=1, keepdims=True) + 1e-9)
    print(f"expression matrix: {len(genes)} genes x {len(sample_cols)} samples")

    # ── sample labels ──────────────────────────────────────────────────
    meta = {}
    for row in csv.DictReader(open(args.samples)):
        meta[row["id.sample name"]] = {
            "flight": row["study.factor value.spaceflight"] == "Space Flight",
            "pd": "Parkinson" in row["study.factor value.donor medical history"],
            "microglia": row["study.factor value.co-culture"] == "with microglia",
        }
    labels = [meta[s] for s in sample_cols]

    # ── embeddings: PPR with up/down entry weights, concatenated ───────
    emb = np.zeros((len(sample_cols), 2 * N))
    for j in range(len(sample_cols)):
        up, dn = np.zeros(N), np.zeros(N)
        for gi, g in enumerate(genes):
            z = Z[gi, j]
            (up if z > 0 else dn)[gene_node[g]] = abs(z)
        for half, vec in ((0, up), (1, dn)):
            if vec.sum() > 0:
                emb[j, half * N:(half + 1) * N] = personalized_pagerank(P, vec)
    np.save(out / "embeddings.npy", emb)
    print(f"embeddings: {emb.shape}")

    # ── leave-one-out classification ───────────────────────────────────
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline

    results = {}
    for task in ("flight", "pd", "microglia"):
        y = np.array([int(l[task]) for l in labels])
        preds = []
        for tr, te in LeaveOneOut().split(emb):
            clf = make_pipeline(StandardScaler(), PCA(n_components=10),
                                LogisticRegression(max_iter=2000, C=1.0))
            clf.fit(emb[tr], y[tr])
            preds.append(clf.predict(emb[te])[0])
        acc = float((np.array(preds) == y).mean())
        base = float(max(y.mean(), 1 - y.mean()))
        results[task] = {"loo_accuracy": round(acc, 3), "majority_baseline": round(base, 3),
                         "n_pos": int(y.sum()), "n": len(y)}
        print(f"  {task:10s} LOO accuracy {acc:.3f} (baseline {base:.3f})")

    # ── PCA coordinates for visualization ──────────────────────────────
    coords = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(emb))
    with (out / "pca_coords.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "pc1", "pc2", "flight", "pd", "microglia"])
        for j, s in enumerate(sample_cols):
            w.writerow([s, round(float(coords[j, 0]), 4), round(float(coords[j, 1]), 4),
                        labels[j]["flight"], labels[j]["pd"], labels[j]["microglia"]])
    json.dump({"n_nodes": N, "n_edges": len(edge_set), "n_genes": len(genes),
               "alpha": ALPHA, "results": results},
              open(out / "results.json", "w"), indent=1)

    # scatter plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        for ax, task, title in zip(axes, ("flight", "pd", "microglia"),
                                   ("Space Flight vs Ground", "PD vs healthy donor",
                                    "with vs without microglia")):
            m = np.array([l[task] for l in labels])
            ax.scatter(coords[~m, 0], coords[~m, 1], c="#2a78d6", s=60, label="no")
            ax.scatter(coords[m, 0], coords[m, 1], c="#eb6834", s=60, label="yes")
            ax.set_title(f"{title}\nLOO acc {results[task]['loo_accuracy']:.2f} "
                         f"(base {results[task]['majority_baseline']:.2f})", fontsize=10)
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=8)
        fig.suptitle("OSD-871 samples in SPOKE embedding space (PSEV-style PPR)", fontsize=12)
        fig.tight_layout()
        fig.savefig(out / "embedding_pca.png", dpi=150)
        print(f"wrote {out}/embedding_pca.png")
    except ImportError:
        print("matplotlib unavailable — skipped plot")


if __name__ == "__main__":
    main()

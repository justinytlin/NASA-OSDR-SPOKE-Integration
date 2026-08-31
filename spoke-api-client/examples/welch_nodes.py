#!/usr/bin/env python3
"""Which SPOKE nodes separate two sample groups?

The meta-analysis idea of Nelson et al. 2021 (Methods 2.5 / replication
welch_meta.py) applied to the per-sample PPR embeddings produced by
spoke_embeddings_ml.py: Welch's t-test per subgraph node between the two
classes of a label, run separately on the up- and down-regulation halves of
the embedding, then per-node-type p-value ranks and top-2.5% flags.

Where the classifier says the classes separate, this says WHICH diseases,
pathways, compounds etc. drive the separation.

Usage:
  welch_nodes.py --dir osd871_ml [--labels pca_coords.csv]
      [--task flight pd microglia] [--top-pct 2.5]

Expects in --dir: subgraph.jsonl, embeddings.npy, and a labels CSV with a
'sample' column plus one boolean column per task (pca_coords.csv works).
Writes welch_<task>.csv into --dir.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_node_meta(subgraph_path):
    """Same insertion order as spoke_embeddings_ml.py builds its node index."""
    node_index, node_meta = {}, []
    for line in open(subgraph_path):
        rec = json.loads(line)
        for n in rec["nodes"]:
            if n["id"] not in node_index:
                node_index[n["id"]] = len(node_meta)
                node_meta.append(n)
    return node_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--labels", default="pca_coords.csv",
                    help="CSV in --dir with 'sample' + boolean task columns")
    ap.add_argument("--task", nargs="*", default=None,
                    help="label columns to test (default: all boolean columns)")
    ap.add_argument("--top-pct", type=float, default=2.5)
    args = ap.parse_args()
    d = Path(args.dir)

    node_meta = load_node_meta(d / "subgraph.jsonl")
    n = len(node_meta)
    emb = np.load(d / "embeddings.npy")
    if emb.shape[1] != 2 * n:
        raise SystemExit(f"embeddings ({emb.shape[1]}) != 2 x subgraph nodes ({n}); "
                         "subgraph.jsonl and embeddings.npy are out of sync")

    labels = pd.read_csv(d / args.labels)
    tasks = args.task or [c for c in labels.columns
                          if set(labels[c].astype(str)) <= {"True", "False"}]
    print(f"{n} nodes, {emb.shape[0]} samples; tasks: {', '.join(tasks)}")

    types = np.array([m["type"] for m in node_meta])
    base = pd.DataFrame({
        "node_id": [m["id"] for m in node_meta],
        "node_type": types,
        "name": [m["name"] for m in node_meta],
    })

    for task in tasks:
        y = labels[task].astype(str).values == "True"
        if y.sum() < 2 or (~y).sum() < 2:
            print(f"[{task}] skipped: needs >=2 samples per class")
            continue
        res = base.copy()
        for half, sl in (("up", slice(0, n)), ("down", slice(n, 2 * n))):
            t, p = stats.ttest_ind(emb[y, sl], emb[~y, sl], axis=0, equal_var=False)
            res[f"t_{half}"] = t
            res[f"p_{half}"] = np.where(np.isnan(p), 1.0, p)
        res["best_p"] = res[["p_up", "p_down"]].min(axis=1)
        res["rank_by_type"] = res.groupby("node_type")["best_p"].rank(method="first")
        thresh = res.groupby("node_type")["node_id"].transform("size") * args.top_pct / 100.0
        res[f"top_{args.top_pct}pct"] = res["rank_by_type"] <= np.ceil(thresh)

        out = d / f"welch_{task}.csv"
        res.sort_values("best_p").to_csv(out, index=False)
        n_top = int(res[f"top_{args.top_pct}pct"].sum())
        print(f"\n[{task}] {y.sum()} vs {(~y).sum()} samples; "
              f"{n_top} top-{args.top_pct}% nodes -> {out.name}")
        for node_type in ("Disease", "Pathway", "BiologicalProcess", "Compound",
                          "Symptom"):
            sub = res[res.node_type == node_type].nsmallest(5, "best_p")
            if not len(sub):
                continue
            print(f"  {node_type}:")
            for _, r in sub.iterrows():
                direction = "up" if r.p_up <= r.p_down else "down"
                print(f"    p={r.best_p:.2e} ({direction})  {r['name']}")


if __name__ == "__main__":
    main()

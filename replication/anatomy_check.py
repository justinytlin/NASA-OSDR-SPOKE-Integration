#!/usr/bin/env python3
"""Tissue sanity check: do a rank table's extreme Anatomy nodes make sense?

The source tissue of a dataset is known ground truth, so Anatomy ranks are a
built-in positive control for the PSEV projection. Because ranks reflect the
DIFFERENTIAL signature (not bulk tissue expression), the source tissue can
legitimately surface at either tail: near rank 1 if its characteristic genes
are up-regulated by the condition, near the bottom if they are suppressed
(e.g. lymphoid organs under spaceflight immune suppression).

Prints the top/bottom Anatomy nodes for each Space-v-Ground-style comparison
(or --col / --col-regex to choose) and, with --tissue, the source tissue
node's mean rank-by-type percentile.

Usage:
  anatomy_check.py RANKS.tsv [--tissue "spleen" ...] [--n 8]
      [--col-regex "SPACE.*GROUND"] [--include-stages]
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd

# UBERON life-stage / fluid / process nodes with sparse gene annotations sit
# at the extremes for spurious reasons; exclude them from the top/bottom view
STAGE_RE = re.compile(
    r"stage|life cycle|zygote|larval|pupal|neonate|embryo|fluid|secretion|"
    r"milk|colostrum|processual|temporal boundary", re.IGNORECASE)


def pick_cols(df, col, col_regex):
    fc_cols = [c for c in df.columns if c.startswith("Log2fc_(")]
    if col:
        return [c for c in fc_cols if c == col or c[len("Log2fc_"):] == col]
    if col_regex:
        return [c for c in fc_cols if re.search(col_regex, c, re.IGNORECASE)]
    out = []
    for c in fc_cols:  # default: Space-v-Ground style comparisons
        p = c[len("Log2fc_("):-1].split(")v(")
        if len(p) == 2 and re.search("SPACE|FLIGHT|FLT", p[0], re.I) \
                and re.search("GROUND|GC", p[1], re.I):
            out.append(c)
    return out or fc_cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ranks", help="*_ranks_and_rank_by_type_*.tsv")
    ap.add_argument("--tissue", action="append", default=[],
                    help="source tissue node name(s) to locate; repeatable")
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--col", default=None, help="exact comparison column")
    ap.add_argument("--col-regex", default=None)
    ap.add_argument("--include-stages", action="store_true",
                    help="keep life-stage/fluid nodes in the top/bottom view")
    args = ap.parse_args()

    df = pd.read_csv(args.ranks, sep="\t", low_memory=False)
    cols = pick_cols(df, args.col, args.col_regex)
    if not cols:
        sys.exit("no comparison columns matched")
    ana = df[df.Node_Type == "Anatomy"].copy()
    ana["mean_bt"] = ana[["Rank_by_type_" + c for c in cols]].mean(axis=1)
    n = len(ana)
    print(f"{len(cols)} comparison(s) pooled; {n} anatomy nodes")

    view = ana if args.include_stages else \
        ana[~ana.Node_Name.astype(str).str.contains(STAGE_RE)]
    print(f"\ntop {args.n} (associated with UP-regulated genes):")
    for _, r in view.nsmallest(args.n, "mean_bt").iterrows():
        print(f"  {r.mean_bt:7.0f} ({100*r.mean_bt/n:5.1f}%)  {r.Node_Name}")
    print(f"\nbottom {args.n} (associated with DOWN-regulated genes):")
    for _, r in view.nlargest(args.n, "mean_bt").iterrows():
        print(f"  {r.mean_bt:7.0f} ({100*r.mean_bt/n:5.1f}%)  {r.Node_Name}")

    for t in args.tissue:
        hit = ana[ana.Node_Name.astype(str).str.lower() == t.lower()]
        if not len(hit):
            near = ana[ana.Node_Name.astype(str).str.contains(t, case=False)]
            print(f"\nsource tissue '{t}': no exact node"
                  + (f"; close names: "
                     + ", ".join(near.Node_Name.astype(str).head(5)) if len(near) else ""))
            continue
        r = hit.iloc[0]
        pct = 100 * r.mean_bt / n
        tail = "UP-associated" if pct < 50 else "DOWN-associated"
        print(f"\nsource tissue '{t}': mean rank {r.mean_bt:.0f}/{n} "
              f"= {pct:.1f}th percentile ({tail} tail; extremes on either "
              f"side = tissue-identity genes shifted by the condition)")


if __name__ == "__main__":
    main()

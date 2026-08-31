#!/usr/bin/env python3
"""Negative control: PSEV ranks from within-group sample splits (no biology).

Takes the per-sample count columns of an OSDR/GeneLab DE table, restricts to
one condition group (e.g. ground control), splits it into random halves,
computes a fake log2FC per gene for each split, and runs every split through
the same PSEV machinery as the real analysis — all splits share one two-pass
sweep of gene_psev.zip, so the whole null ensemble costs one run.

Nodes that rank high in these no-biology nulls (well-connected hubs) rank
high for any large gene set; a real finding should beat the null ensemble.
With --real REAL_RANKS.tsv the script writes a per-node comparison with an
empirical percentile of the real rank against the null splits.

Usage:
  osdr_null_control.py TABLE.csv --name OSD-564-GCnull --sample-pattern _GC_
      [--human] [--n-splits 10] [--seed 0]
      [--n-genes 1415 | --p 0.01]       per-split gene selection (top-N by
                                        Welch p between halves, or cutoff)
      [--real RANKS.tsv --real-col "Log2fc_(...)v(...)"]
      [--out DIR] [--zip ...] [--spoke ...] [--homologene ...] [-b 0.1]
      [--prep-only] [--limit-groups N]
"""

import argparse
import itertools
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(HERE))
from osdr_to_psev import load_de_table, map_to_spoke_genes
from psev_pipeline import compute_rank_table


def make_splits(n, n_splits, rng):
    """Half-splits of range(n) as (A, B) index tuples, complements deduped."""
    k = n // 2
    combos = [c for c in itertools.combinations(range(n), k)
              if n % 2 == 1 or 0 in c]          # even n: fix sample 0 in A
    if len(combos) > n_splits:
        combos = [combos[i] for i in rng.choice(len(combos), n_splits, replace=False)]
    return [(list(c), [i for i in range(n) if i not in c]) for c in combos]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--name", required=True)
    ap.add_argument("--sample-pattern", required=True,
                    help="regex selecting the control group's sample columns, e.g. _GC_")
    ap.add_argument("--human", action="store_true")
    ap.add_argument("--n-splits", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-genes", type=int, default=None,
                    help="per split keep the N genes with smallest Welch p "
                         "(match the real run's signature size)")
    ap.add_argument("--p", type=float, default=0.01,
                    help="p cutoff per split when --n-genes is not given")
    ap.add_argument("--real", default=None,
                    help="real rank table to compare against (osdr_to_psev output)")
    ap.add_argument("--real-col", default=None,
                    help="fc column of the real run, e.g. 'Log2fc_(Space Flight)v(Ground Control)'")
    ap.add_argument("--out", default=str(HERE / "psev_out"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    ap.add_argument("--prep-only", action="store_true")
    ap.add_argument("--limit-groups", type=int, default=None)
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)
    save_str = "_".join(str(args.b).split("."))

    node_info_df = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    spoke_genes = set(node_info_df[node_info_df.Node_Type == "Gene"].Node.values)

    df = load_de_table(args.table)
    sample_cols = [c for c in df.columns
                   if pd.Series([c]).str.contains(args.sample_pattern, regex=True)[0]
                   and pd.api.types.is_numeric_dtype(df[c])]
    if len(sample_cols) < 4:
        sys.exit(f"only {len(sample_cols)} sample columns match "
                 f"'{args.sample_pattern}' — need at least 4 to split")
    print(f"[{args.name}] {len(sample_cols)} control samples: "
          f"{', '.join(sample_cols)}")

    # the per-sample columns in OSDR DE tables are UNnormalized counts;
    # library sizes can differ >2x, so apply DESeq2-style median-of-ratios
    # size factors before computing split fold changes
    sub = df[["ENTREZID"] + sample_cols].copy()
    raw = sub[sample_cols].astype(float).values
    pos = (raw > 0).all(axis=1)
    lg = np.log(raw[pos])
    sf = np.exp(np.median(lg - lg.mean(axis=1)[:, None], axis=0))
    print("  size factors: " + ", ".join(f"{s:.2f}" for s in sf))
    sub[sample_cols] = np.log2(raw / sf + 1.0)
    mapped = map_to_spoke_genes(sub, spoke_genes, args.homologene, human=args.human)
    expr = mapped.groupby("Node")[sample_cols].mean()
    X = expr.values
    print(f"  mapped SPOKE genes: {len(expr)}")

    rng = np.random.default_rng(args.seed)
    splits = make_splits(len(sample_cols), args.n_splits, rng)
    print(f"  within-group splits: {len(splits)} "
          f"({len(splits[0][0])} vs {len(splits[0][1])} samples)")

    kept = pd.DataFrame({"Node": expr.index})
    split_cols = []
    for si, (a, b) in enumerate(splits):
        fc = X[:, a].mean(axis=1) - X[:, b].mean(axis=1)
        with warnings.catch_warnings():
            # constant-count genes trigger scipy precision warnings; their p
            # comes back nan and is set to 1 below
            warnings.simplefilter("ignore", RuntimeWarning)
            _, p = stats.ttest_ind(X[:, a], X[:, b], axis=1, equal_var=False)
        p = np.where(np.isnan(p), 1.0, p)
        if args.n_genes:
            sel = np.argsort(p)[:args.n_genes]
            keep = np.zeros(len(fc), dtype=bool)
            keep[sel] = True
        else:
            keep = p < args.p
        col = f"Nullfc_split{si:02d}"
        kept[col] = np.where(keep, fc, 0.0)
        split_cols.append(col)
        print(f"  split {si:02d}: {int(keep.sum())} genes selected")
    kept = kept[np.abs(kept[split_cols].values).max(axis=1) > 0]
    print(f"  union of selected genes: {len(kept)}")
    if args.prep_only:
        return

    result_df = compute_rank_table(kept, split_cols, node_info_df, args.zip,
                                   b=args.b, limit_groups=args.limit_groups)
    out = out_dir / f"{args.name}_null_ranks_and_rank_by_type_for_meta{save_str}.tsv"
    result_df.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out}")

    if not args.real:
        return
    if not args.real_col:
        sys.exit("--real given without --real-col")
    real = pd.read_csv(args.real, sep="\t")
    bt_real = "Rank_by_type_" + args.real_col
    null_bt = result_df[[f"Rank_by_type_{c}" for c in split_cols]].values
    S = null_bt.shape[1]
    cmp_df = node_info_df.copy()
    cmp_df["real_rank_by_type"] = real[bt_real].values
    cmp_df["null_median"] = np.median(null_bt, axis=1)
    cmp_df["null_best"] = null_bt.min(axis=1)
    # empirical percentile of the real rank vs the null ensemble, both tails
    cmp_df["p_top"] = (1 + (null_bt <= cmp_df.real_rank_by_type.values[:, None]).sum(axis=1)) / (1 + S)
    cmp_df["p_bottom"] = (1 + (null_bt >= cmp_df.real_rank_by_type.values[:, None]).sum(axis=1)) / (1 + S)
    type_size = cmp_df.groupby("Node_Type")["Node"].transform("size")
    cmp_df["real_top2.5pct"] = cmp_df.real_rank_by_type <= np.ceil(type_size * 0.025)
    cmp_df["beats_all_nulls"] = cmp_df["real_top2.5pct"] & (cmp_df.p_top == 1 / (1 + S))
    cmp_out = out_dir / f"{args.name}_vs_real.tsv"
    cmp_df.to_csv(cmp_out, sep="\t", index=False)
    print(f"  wrote {cmp_out}")

    top = cmp_df[cmp_df["real_top2.5pct"]]
    print(f"\n  real top-2.5% nodes: {len(top)}; "
          f"beating every null split: {int(top.beats_all_nulls.sum())}")
    for t in ("Disease", "Pathway", "BiologicalProcess", "Symptom", "Anatomy"):
        sub = top[top.Node_Type == t]
        if not len(sub):
            continue
        n_rob = int(sub.beats_all_nulls.sum())
        print(f"  {t}: {n_rob}/{len(sub)} robust; flagged as null-prone e.g.: "
              + "; ".join(sub[~sub.beats_all_nulls].nsmallest(3, "real_rank_by_type")
                          .Node_Name.astype(str)))


if __name__ == "__main__":
    main()

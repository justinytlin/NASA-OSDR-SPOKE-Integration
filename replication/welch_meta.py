#!/usr/bin/env python3
"""Meta-analysis step of Nelson et al., Life 2021 (Methods 2.5).

Pools per-comparison PSEV node ranks from the six GLDS studies into three
groups (Ground v Basal, Space v Basal, Space v Ground), runs Welch's t-test
per node (space groups vs Ground-v-Basal), selects the top 2.5% per node type,
and compares the outcome with the paper's shipped ground truth
(top_welch_nodes.xlsx).

Usage:
  welch_meta.py [--rank-col overall|bytype] [--direction num|den] [--out DIR]
"""

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
V2 = Path(os.environ.get("WELCH_V2", BASE / "SPOKE_NASA/GeneLab_for_SPOKE/V2"))
ACCESSIONS = ["GLDS-4", "GLDS-244", "GLDS-245", "GLDS-246", "GLDS-288", "GLDS-289"]

GROUP_RE = re.compile(r"SPACE|GROUND|BASAL", re.IGNORECASE)


def classify(comparison_col):
    """'Log2fc_(A)v(B)' -> ('Space','Ground') style tuple or None."""
    inner = comparison_col[len("Log2fc_("):-1]
    parts = inner.split(")v(")
    if len(parts) != 2:
        return None
    labels = []
    for p in parts:
        m = GROUP_RE.search(p)
        if not m:
            return None
        labels.append(m.group(0).capitalize())
    return tuple(labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank-col", choices=["overall", "bytype"], default="overall")
    ap.add_argument("--direction", choices=["num", "den", "both"], default="num",
                    help="num: use (Space)v(Ground) columns for SvG; den: use (Ground)v(Space)")
    ap.add_argument("--ag", choices=["space", "ground", "drop"], default="space",
                    help="how to treat 'Space Flight & 1G by centrifugation' groups")
    ap.add_argument("--matched", action="store_true",
                    help="only pool comparisons with matched mission/timepoint/collection")
    ap.add_argument("--out", default=str(HERE / "meta_out"))
    ap.add_argument("--suffix", default="0_1")
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"),
                    help="dir holding node_info.tsv of the PSEV universe the rank tables were built on")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(exist_ok=True)

    def group_of(pair, col=""):
        if len(pair) != 2:
            return None
        # optionally reinterpret artificial-gravity groups
        if "centrifugation" in col:
            inner = col[len("Log2fc_("):-1]
            parts = inner.split(")v(")
            labs = list(pair)
            for i, p in enumerate(parts):
                if "centrifugation" in p:
                    if args.ag == "drop":
                        return None
                    if args.ag == "ground":
                        labs[i] = "Ground"
            pair = tuple(labs)
        a, b = pair
        if a == b:
            return None
        for aa, bb in ((a, b), (b, a)) if args.direction == "both" else (((b, a),) if args.direction == "den" else ((a, b),)):
            if (aa, bb) == ("Space", "Ground"):
                return "SvG"
            if (aa, bb) == ("Space", "Basal"):
                return "SvB"
            if (aa, bb) == ("Ground", "Basal"):
                return "GvB"
        return None

    dur_re = re.compile(r"^~?\d+( day)?$")

    def is_matched(col):
        inner = col[len("Log2fc_("):-1]
        a, b = inner.split(")v(")
        fa, fb = a.split(" & "), b.split(" & ")
        # mission (MHU-x) must agree
        ma = [f for f in fa if f.startswith("MHU")]
        mb = [f for f in fb if f.startswith("MHU")]
        if ma != mb:
            return False
        # collection protocol must agree when both sides carry one
        coll = ("Carcass", "Upon euthanasia")
        ca = [f for f in fa if f in coll]
        cb = [f for f in fb if f in coll]
        if ca and cb and ca != cb:
            return False
        # duration must agree unless one side is the basal group
        if "BASAL" not in a.upper() and "BASAL" not in b.upper():
            da = [f for f in fa if dur_re.match(f)]
            db = [f for f in fb if dur_re.match(f)]
            if da != db:
                return False
        return True

    node_info = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t", dtype={"Node": str})
    groups = {"GvB": [], "SvB": [], "SvG": []}

    for acc in ACCESSIONS:
        path = V2 / f"{acc}_ranks_and_rank_by_type_for_meta{args.suffix}.tsv"
        if not path.exists():
            print(f"!! missing {path.name}, skipping")
            continue
        df = pd.read_csv(path, sep="\t")
        fc_cols = [c for c in df.columns if c.startswith("Log2fc_(")]
        for col in fc_cols:
            g = group_of(classify(col) or (), col)
            if g is None or (args.matched and not is_matched(col)):
                continue
            use = col if args.rank_col == "overall" else f"Rank_by_type_{col}"
            groups[g].append(df[use].values.astype(np.float64))
        print(f"{acc}: " + ", ".join(
            f"{g}={sum(1 for c in fc_cols if group_of(classify(c) or (), c) == g)}"
            for g in groups))

    for g, cols in groups.items():
        print(f"group {g}: {len(cols)} comparisons")
    gvb = np.array(groups["GvB"])
    result = node_info.copy()

    for label, key in (("Space v Basal", "SvB"), ("Space v Ground", "SvG")):
        mat = np.array(groups[key])
        t, p = stats.ttest_ind(mat, gvb, axis=0, equal_var=False)
        result[f"P-value {label} - Ground v Basal"] = p
        result[f"T-stat {label} - Ground v Basal"] = t
    both = np.vstack([groups["SvB"], groups["SvG"]])
    t, p = stats.ttest_ind(both, gvb, axis=0, equal_var=False)
    result["P-value combined space"] = p

    # per-type p-value ranks and top 2.5% flags
    for label in ("Space v Basal", "Space v Ground"):
        col = f"P-value {label} - Ground v Basal"
        result[f"Rank by type {label}"] = (
            result.groupby("Node_Type")[col].rank(method="first"))
        thresh = result.groupby("Node_Type")["Node"].transform("size") * 0.025
        result[f"Top 2.5% {label}"] = result[f"Rank by type {label}"] <= np.ceil(thresh)
    result["Top 2.5% either space comparison"] = (
        result["Top 2.5% Space v Basal"] | result["Top 2.5% Space v Ground"])

    n_top = int(result["Top 2.5% either space comparison"].sum())
    print(f"top nodes (2.5% either comparison): {n_top} "
          f"({100 * n_top / len(result):.1f}%)  [paper: 15,801; 4.1%]")

    result.to_csv(out / f"welch_results_{args.rank_col}_{args.direction}.tsv",
                  sep="\t", index=False)

    # ---- compare with paper ground truth --------------------------------
    truth = pd.read_excel(BASE / "SPOKE_NASA/top_welch_nodes.xlsx",
                          sheet_name="all top nodes")
    top25_col = [c for c in truth.columns if c.startswith("Top 2.5%")]
    if top25_col:
        strict = truth[truth[top25_col[0]] == True]
        if len(strict):
            truth = strict
    truth_nodes = set(truth.Node.astype(str))
    ours = set(result[result["Top 2.5% either space comparison"]].Node.astype(str))
    inter = truth_nodes & ours
    print(f"paper top nodes: {len(truth_nodes)}; ours: {len(ours)}; "
          f"overlap: {len(inter)} ({100 * len(inter) / max(len(truth_nodes), 1):.1f}% of paper's)")

    merged = truth.assign(Node=truth.Node.astype(str)).merge(
        result.assign(Node=result.Node.astype(str)),
        on="Node", suffixes=(" (paper)", " (ours)"))
    for paper_col, our_col in [
        ("P-value Space v Ground - Ground v Basal", "P-value Space v Ground - Ground v Basal (ours)"),
        ("P-value Space v Basal - Ground v Basal", "P-value Space v Basal - Ground v Basal (ours)"),
    ]:
        pc = paper_col + " (paper)" if paper_col + " (paper)" in merged.columns else paper_col
        if pc in merged.columns and our_col in merged.columns:
            rho = stats.spearmanr(merged[pc], merged[our_col]).statistic
            print(f"Spearman rho ({paper_col}): {rho:.3f}")

    # the paper's 22 showcased nodes (Figure 4)
    showcase = ["Space Motion Sickness", "taste receptor complex",
                "Vitamin D (calciferol) metabolism", "sympathetic nervous system",
                "regulation of vasoconstriction", "regulation of blood vessel diameter",
                "regulation of cortisol secretion"]
    for name in showcase:
        row = result[result.Node_Name.astype(str).str.lower() == name.lower()]
        if len(row):
            r = row.iloc[0]
            print(f"  {name}: top2.5%={bool(r['Top 2.5% either space comparison'])} "
                  f"pSvG={r['P-value Space v Ground - Ground v Basal']:.4g} "
                  f"pSvB={r['P-value Space v Basal - Ground v Basal']:.4g}")


if __name__ == "__main__":
    main()

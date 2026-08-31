#!/usr/bin/env python3
"""Reproduce Figure 4 violin panels: node-rank distributions per comparison group.

For each showcase node, collect its overall PSEV rank across the pooled FC
comparisons in each group (Ground v Basal, Space v Basal, Space v Ground) and
draw violins, matching the paper's colors (blue/yellow/green).

Usage: plot_violins.py [--out FILE.png]
"""

import argparse
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
V2 = Path(os.environ.get("WELCH_V2", BASE / "SPOKE_NASA/GeneLab_for_SPOKE/V2"))
ACCESSIONS = ["GLDS-4", "GLDS-244", "GLDS-245", "GLDS-246", "GLDS-288", "GLDS-289"]
GROUP_RE = re.compile(r"SPACE|GROUND|BASAL", re.IGNORECASE)

SHOWCASE = [
    # (Node_Name as in node_info, Node_Type)
    ("Polydipsia, Psychogenic", "Symptom"),
    ("Space Motion Sickness", "Symptom"),
    ("Hearing Loss, Conductive", "Symptom"),
    ("Vision, Low", "Symptom"),
    ("Deaf-Blind Disorders", "Symptom"),
    ("Tachypnea", "Symptom"),
    ("taste receptor complex", "CellularComponent"),
    ("regulation of cortisol secretion", "BiologicalProcess"),
    ("Vitamin D (calciferol) metabolism", "Pathway"),
    ("regulation of vasoconstriction", "BiologicalProcess"),
    ("regulation of blood vessel diameter", "BiologicalProcess"),
    ("sympathetic nervous system", "Anatomy"),
    ("autonomic nervous system", "Anatomy"),
    ("corticobulbar and corticospinal tracts", "Anatomy"),
    ("pyramidal decussation", "Anatomy"),
    ("hindbrain arachnoid mater", "Anatomy"),
    ("spinal cord ependyma", "Anatomy"),
]
COLORS = {"GvB": "#2b6ca3", "SvB": "#f2a534", "SvG": "#3d9142"}
LABELS = {"GvB": "Ground v Baseline", "SvB": "Space v Baseline", "SvG": "Space v Ground"}


def classify(col):
    inner = col[len("Log2fc_("):-1]
    parts = inner.split(")v(")
    if len(parts) != 2:
        return None
    labs = []
    for p in parts:
        m = GROUP_RE.search(p)
        if not m:
            return None
        labs.append(m.group(0).capitalize())
    a, b = labs
    return {"Space,Ground": "SvG", "Space,Basal": "SvB",
            "Ground,Basal": "GvB"}.get(f"{a},{b}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "meta_out/figure4_violins.png"))
    ap.add_argument("--suffix", default="0_1")
    args = ap.parse_args()

    ranks = {g: [] for g in ("GvB", "SvB", "SvG")}
    node_meta = None
    for acc in ACCESSIONS:
        path = V2 / f"{acc}_ranks_and_rank_by_type_for_meta{args.suffix}.tsv"
        if not path.exists():
            print(f"!! missing {path.name}")
            continue
        df = pd.read_csv(path, sep="\t")
        if node_meta is None:
            node_meta = df[["Node", "Node_Name", "Node_Type"]]
        for col in df.columns:
            if col.startswith("Log2fc_("):
                g = classify(col)
                if g:
                    ranks[g].append(df[col].values.astype(float))
    mats = {g: np.array(v) for g, v in ranks.items()}

    present = []
    for name, ntype in SHOWCASE:
        sel = (node_meta.Node_Name.astype(str) == name) & (node_meta.Node_Type == ntype)
        if sel.sum():
            present.append((name, ntype, np.flatnonzero(sel.values)[0]))
        else:
            print(f"!! node not found: {name} ({ntype})")

    ncols = 6
    nrows = int(np.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows))
    for ax in np.ravel(axes):
        ax.axis("off")
    for k, (name, ntype, i) in enumerate(present):
        ax = np.ravel(axes)[k]
        ax.axis("on")
        data = [mats[g][:, i] for g in ("GvB", "SvB", "SvG")]
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for body, g in zip(parts["bodies"], ("GvB", "SvB", "SvG")):
            body.set_facecolor(COLORS[g])
            body.set_alpha(0.8)
        for elem in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[elem].set_color("#333333")
        ax.set_title(f"{name}\n({ntype})", fontsize=9)
        ax.set_xticks([])
        ax.set_ylabel("rank", fontsize=8)
        ax.tick_params(labelsize=7)
    handles = [plt.Line2D([0], [0], color=COLORS[g], lw=6, label=LABELS[g])
               for g in ("GvB", "SvB", "SvG")]
    fig.legend(handles=handles, loc="lower right", fontsize=10)
    fig.suptitle("Replication of Nelson et al. 2021 Fig. 4 violin panels "
                 "(pooled PSEV node ranks by comparison group)", fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    Path(args.out).parent.mkdir(exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

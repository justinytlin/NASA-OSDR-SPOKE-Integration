#!/usr/bin/env python3
"""Compile the per-dataset results of a build run (e.g. the 2025 upregulates-only build).

Reads, for every entry of data/_catalog/manifest_<tag>.json, datasets/<name>/corrected_linked_<tag>/summary.json
(rank of the tissue's own gene-linked node) and datasets/<name>/corrected_ontology/summary.json (ontology-family
check), and writes results/compiled_<tag>.tsv, results/tissue_summary_<tag>.tsv and
results/tissue_overview_<tag>.png (one row per tissue: own-node rank per dataset, verdict-coded).
Usage: compile_build.py --tag 2025upAuG
"""
import argparse, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; sys.path.insert(0, str(HERE)); import style as S
VERDICT = {"CLEAR": ("#0ca30c", "o"), "PARTIAL": ("#fab219", "^"), "NONE": ("#d03b3b", "X"), "NO LINKED NODE": (S.GRAY_MARK, "s")}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True)
    ap.add_argument("--linked-tag", default=None, help="folder suffix of the linked check to compile (default: --tag)"); a = ap.parse_args()
    ltag = a.linked_tag or a.tag
    rows = []
    for e in json.load(open(ROOT / "data/_catalog" / f"manifest_{a.tag}.json")):
        d = ROOT / "datasets" / e["name"]
        lj = d / f"corrected_linked_{ltag}" / "summary.json"; oj = d / "corrected_ontology" / "summary.json"
        if not lj.exists():
            continue
        L = json.load(open(lj)); O = json.load(open(oj)) if oj.exists() else {}
        rows.append(dict(dataset=e["base"], payload=e["payload"], tissue=e["tissue"], family=e["family"], control_type=e["control_type"],
                         n_samples=L.get("n_samples"), n_linked=L.get("n_linked"), proxy=bool(L.get("proxy")),
                         linked_verdict=L["verdict"], own_node_rank=L.get("best_target_rank"),
                         frac_samples_rank1=L.get("frac_samples_target_rank1"), target_nodes="; ".join(L.get("target_nodes", [])[:2]),
                         top_linked_nodes="; ".join(L.get("top_linked_nodes", [])[:4]),
                         onto_verdict=O.get("verdict"), onto_family_median_pct=round(O["family_median_pct"], 1) if O else None,
                         onto_family_rank=O.get("family_rank_among_organ_families"), onto_mwu_p=O.get("mwu_p_less")))
    df = pd.DataFrame(rows).sort_values(["tissue", "own_node_rank"]); out = ROOT / "results"
    df.to_csv(out / f"compiled_{ltag}.tsv", sep="\t", index=False)
    g = df.groupby("tissue")
    ts = g.agg(n_datasets=("dataset", "size"), target=("target_nodes", lambda v: re.sub(r" \(rank \d+\)", "", str(v.iloc[0])).split(";")[0]),
               proxy=("proxy", "first"), linked_CLEAR=("linked_verdict", lambda v: int((v == "CLEAR").sum())),
               linked_PARTIAL=("linked_verdict", lambda v: int((v == "PARTIAL").sum())), linked_NONE=("linked_verdict", lambda v: int((v == "NONE").sum())),
               no_linked_node=("linked_verdict", lambda v: int((v == "NO LINKED NODE").sum())),
               median_own_rank=("own_node_rank", "median"), onto_CLEAR=("onto_verdict", lambda v: int((v == "CLEAR").sum())),
               onto_PARTIAL=("onto_verdict", lambda v: int((v == "PARTIAL").sum())), onto_NONE=("onto_verdict", lambda v: int((v == "NONE").sum())),
               onto_median_family_pct=("onto_family_median_pct", "median")).reset_index().sort_values("median_own_rank")
    ts.to_csv(out / f"tissue_summary_{ltag}.tsv", sep="\t", index=False)
    # overview figure
    S.setup(); order = ts.tissue.tolist(); n_linked = int(df.n_linked.dropna().iloc[0]) if df.n_linked.notna().any() else 44
    h = 0.33 * len(order) + 2.0
    fig, ax = plt.subplots(figsize=(9, h)); fig.subplots_adjust(left=0.27, right=0.97, top=1 - 0.5 / h, bottom=1.25 / h)
    yy = np.arange(len(order))[::-1]; rng = np.random.default_rng(0)
    for y, t in zip(yy, order):
        d = df[df.tissue == t]
        for v, (c, mk) in VERDICT.items():
            gg = d[d.linked_verdict == v]
            x = gg.own_node_rank.fillna(n_linked + 4).values
            nolink = v == "NO LINKED NODE"
            ax.scatter(x, y + rng.uniform(-0.16, 0.16, len(gg)), s=46 if not nolink else 34, marker=mk, zorder=3,
                       facecolors="white" if nolink else c, edgecolors=S.INK2 if nolink else "white", linewidths=1.0 if nolink else 0.8)
    ax.axvline(n_linked + 4, color=S.GRID, lw=0.8, ls=(0, (3, 3)))
    ax.set_yticks(yy); ax.set_yticklabels([f"{t} ({n})" + (" proxy" if p else "") for t, n, p in zip(ts.tissue, ts.n_datasets, ts.proxy)], fontsize=9)
    ax.set_xlim(0, n_linked + 7); ax.set_xticks([1, 10, 20, 30, 40, n_linked + 4]); ax.set_xticklabels(["1", "10", "20", "30", "40", "no linked\nnode"], fontsize=9)
    ax.set_xlabel(f"rank of the tissue's own gene-linked Anatomy node among {n_linked} linked nodes (1 = top)", fontsize=9)
    ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], marker=mk, ls="", color=c, mfc="white" if v == "NO LINKED NODE" else c, mec=S.INK2 if v == "NO LINKED NODE" else c, ms=7, label=v.lower()) for v, (c, mk) in VERDICT.items()],
              loc="upper center", bbox_to_anchor=(0.5, -1.0 / (0.33 * len(order))), ncol=4, frameon=False, fontsize=9)
    ax.set_title(f"Own-tissue node rank, {ltag}: {len(df)} datasets, {len(order)} tissues", fontsize=12, pad=10)
    fig.savefig(out / f"tissue_overview_{ltag}.png", dpi=180); fig.savefig(out / f"tissue_overview_{ltag}.svg")
    pd.set_option("display.width", 250)
    print(f"{len(df)} datasets;", df.linked_verdict.value_counts().to_dict(), "| ontology:", df.onto_verdict.value_counts().to_dict())
    print(ts.to_string(index=False, max_colwidth=28))


if __name__ == "__main__":
    main()

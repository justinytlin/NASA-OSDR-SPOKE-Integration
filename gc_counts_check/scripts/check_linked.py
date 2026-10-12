#!/usr/bin/env python3
"""Tissue identification among the Anatomy nodes that carry gene-expression links.

In the 2019 SPOKE only 74 of the 13,257 Anatomy nodes are linked to genes (found from
the PSEVs themselves by extract_anatomy_psev.py: a linked node has a PSEV value above
1e-4 for more than half of all genes; the next node down has ~6 %). Every other Anatomy
node gets its score second-hand through the ontology. This check asks the direct
question: among those 74 linked tissues, where does the sample's own tissue rank?

Target nodes for a tissue = the linked nodes inside its UBERON family (root term and
everything below it). If the family has none (spleen, thymus, skeletal muscle, colon ...),
the nearest linked ancestors of the root are used instead and the result is flagged
"proxy" — the tissue itself is not represented in the graph's expression data.

Scores are the baseline-corrected scores (scores - mean_weight * uniform); ranks are
within the 74 linked nodes, per sample, then averaged (pooled).
Verdict: CLEAR = a target node is ranked 1st of 74 on the pooled ranking;
PARTIAL = best target in the top 7 (top 10 %); NONE otherwise. Proxy results carry
the suffix "(proxy)".

Output datasets/<name>/corrected_linked/: A_linked_nodes.png/.svg (the 74 nodes ranked,
targets highlighted), C_per_sample.png/.svg (rank of the best target in each sample),
summary.json, linked_nodes.tsv.

Usage: check_linked.py --npz psev_out/<name>_GC_ranks.npz --tissue spleen [--label spleen] --name OSD-246
"""
import argparse, json, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import style as S
from ontology_families import build_families, load_obo, ROOTS


def linked_mask(path):
    z = np.load(path)
    frac = (z["values"] > 1e-4).mean(axis=0)
    return frac > 0.5, frac


def nearest_linked_ancestors(root_ids, linked_ids):
    _, parents, _ = load_obo()
    frontier, seen = set(root_ids), set(root_ids)
    while frontier:
        nxt = set()
        for t in frontier:
            nxt |= parents.get(t, set())
        nxt -= seen
        hit = nxt & linked_ids
        if hit:
            return hit
        seen |= nxt; frontier = nxt
    return set()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True); ap.add_argument("--tissue", required=True); ap.add_argument("--label", default=None)
    ap.add_argument("--name", default=None); ap.add_argument("--out", default=str(ROOT / "datasets"))
    ap.add_argument("--spoke", default=str(ROOT.parent / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--score-key", choices=["scores", "zscore", "perm"], default="scores")
    ap.add_argument("--anatomy-psev", default=str(ROOT / "psev_out/_anatomy_psev.npz"),
                    help="extract_anatomy_psev.py output for the same build as --spoke / the npz")
    ap.add_argument("--tag", default="", help="suffix for the output folder, e.g. 2025 -> corrected_linked_2025")
    ap.add_argument("--linked-ids", default=None,
                    help="TSV with an 'id' column of UBERON ids: use these as the gene-linked set instead of detecting it from the PSEVs "
                         "(e.g. results/aug_linked_nodes.tsv = the 72 Anatomy nodes with UPREGULATES_AuG edges in the live graph)")
    a = ap.parse_args()
    name = a.name or Path(a.npz).stem.split("_")[0]; tissue = a.tissue; label = a.label or tissue
    if tissue not in ROOTS:
        sys.exit(f"no UBERON root for '{tissue}'")
    out = Path(a.out) / name / (("corrected_linked" if a.score_key == "scores" else f"corrected_{a.score_key}_linked") + (f"_{a.tag}" if a.tag else "")); out.mkdir(parents=True, exist_ok=True)
    color = S.TISSUE_COLOR.get(tissue, "#2a78d6")

    ni = pd.read_csv(Path(a.spoke) / "node_info.tsv", sep="\t"); amask = (ni.Node_Type == "Anatomy").values
    an = ni[amask].reset_index(drop=True)
    if a.linked_ids:
        ids = set(pd.read_csv(a.linked_ids, sep="\t").id.astype(str))
        link = an.Node.astype(str).isin(ids).values
    else:
        link, _ = linked_mask(a.anatomy_psev)
    ln = an[link].reset_index(drop=True); L = len(ln)
    z = np.load(a.npz)
    if a.score_key == "zscore":
        corr = z["scores_z"].astype(np.float64) - z["mean_weight"][:, None] * z["uniform_z"].astype(np.float64)[None, :]
        cl = corr[:, amask][:, link]
    elif a.score_key == "perm":
        corr = z["scores"].astype(np.float64) - z["mean_weight"][:, None] * z["uniform"].astype(np.float64)[None, :]
        obs = corr[:, amask]
        P = z["perm_scores"].astype(np.float64) - z["mean_weight"][None, :, None] * z["uniform"].astype(np.float64)[amask][None, None, :]
        zz = np.nan_to_num((obs - P.mean(axis=0)) / np.where(P.std(axis=0) > 0, P.std(axis=0), np.nan))
        cl = zz[:, link]
    else:
        corr = z["scores"].astype(np.float64) - z["mean_weight"][:, None] * z["uniform"].astype(np.float64)[None, :]
        cl = corr[:, amask][:, link]                                # samples x linked nodes
    ranks = (-cl).argsort(axis=1).argsort(axis=1) + 1               # 1 = top, per sample
    pooled = ranks.mean(axis=0); order = np.argsort(pooled); pooled_rank = np.empty(L, int); pooled_rank[order] = np.arange(1, L + 1)
    samples = [str(s) for s in z["samples"]]

    fams, roots = build_families(an.Node.values)
    linked_ids = set(ln.Node.values)
    direct = fams[tissue] & linked_ids
    proxy = not direct
    targets = direct or nearest_linked_ancestors(roots[tissue], linked_ids)
    tmask = ln.Node.isin(targets).values
    if not tmask.any():
        summary = dict(name=name, tissue=tissue, tissue_label=label, n_linked=L, verdict="NO LINKED NODE", target_nodes=[], proxy=True)
        (out / "summary.json").write_text(json.dumps(summary, indent=2)); print(f"{name} ({tissue}) [linked]: NO LINKED NODE"); return
    best = int(pooled_rank[tmask].min())
    best_per_sample = ranks[:, tmask].min(axis=1)
    verdict = "CLEAR" if best == 1 else "PARTIAL" if best <= 7 else "NONE"
    tab = pd.DataFrame({"Node": ln.Node, "Node_Name": ln.Node_Name, "pooled_rank": pooled_rank, "mean_rank": pooled.round(2),
                        "target": tmask}).sort_values("pooled_rank")
    tab.to_csv(out / "linked_nodes.tsv", sep="\t", index=False)
    summary = dict(name=name, tissue=tissue, tissue_label=label, n_samples=len(samples), n_linked=L, proxy=bool(proxy),
                   target_nodes=[f"{r.Node_Name} (rank {r.pooled_rank})" for r in tab[tab.target].itertuples()],
                   best_target_rank=best, frac_samples_target_rank1=float((best_per_sample == 1).mean()),
                   frac_samples_target_top3=float((best_per_sample <= 3).mean()), verdict=verdict,
                   top_linked_nodes=[r.Node_Name for r in tab.head(8).itertuples()])
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    S.setup()
    ptxt = " (proxy)" if proxy else ""
    # A. the linked nodes ranked, top 20 plus any target below
    show = tab.head(20)
    extra = tab[tab.target & (tab.pooled_rank > 20)]
    show = pd.concat([show, extra])
    fig, ax = plt.subplots(figsize=(7.6, 0.3 * len(show) + 1.6)); h = 0.3 * len(show) + 1.6
    fig.subplots_adjust(left=0.36, right=0.96, top=1 - 0.55 / h, bottom=0.95 / h)
    yy = np.arange(len(show))[::-1]
    ax.hlines(yy, 0, show.mean_rank, color=S.GRID, lw=1.2, zorder=1)
    ax.scatter(show.mean_rank, yy, s=[52 if t else 26 for t in show.target], color=[color if t else S.GRAY_MARK for t in show.target],
               edgecolors="white", linewidths=1, zorder=3)
    ax.set_yticks(yy); ax.set_yticklabels([f"{r.pooled_rank}. {r.Node_Name}" for r in show.itertuples()], fontsize=8.5)
    for lab, t in zip(ax.get_yticklabels(), show.target):
        if t:
            lab.set_color(S.INK); lab.set_fontweight("semibold")
    ax.set_xlim(0, L + 1); ax.set_xlabel(f"mean rank across samples among the {L} gene-linked Anatomy nodes (1 = top)", fontsize=9)
    ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=color, ms=7, label=f"{label} node" + (" (nearest linked ancestor)" if proxy else "")),
                       Line2D([], [], marker="o", ls="", color=S.GRAY_MARK, ms=5, label="other gene-linked nodes")],
              loc="upper center", bbox_to_anchor=(0.3, -0.62 / (0.3 * len(show))), ncol=2, frameon=False, fontsize=8.5)
    fig.suptitle(f"{name.split('_')[0]} {label}: gene-linked Anatomy nodes ranked{ptxt}", x=0.5, y=1 - 0.2 / h, fontsize=11.5)
    fig.savefig(out / "A_linked_nodes.png", dpi=180); fig.savefig(out / "A_linked_nodes.svg"); plt.close(fig)

    # C. per sample
    n_s = len(samples)
    fig, ax = plt.subplots(figsize=(max(6.4, 0.32 * n_s + 3.2), 5.2)); fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.32)
    ax.scatter(np.arange(n_s), best_per_sample, s=40, color=color, edgecolors="white", linewidths=1, zorder=3)
    ax.axhline(best, color=S.INK2, lw=1, zorder=2)
    ax.set_xticks(np.arange(n_s)); ax.set_xticklabels([re.sub(r"^(Mmus_C57-6[A-Z]*_[A-Z]+_|RR10_[A-Z]+_)", "", s)[:22] for s in samples], rotation=90, fontsize=7)
    ax.set_ylim(L + 1, 0); ax.set_ylabel(f"rank of the best {label} node among {L} linked nodes", fontsize=9)
    ax.grid(axis="y", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax, keep=("bottom", "left"))
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=color, ms=7, label="sample"), Line2D([], [], color=S.INK2, lw=1, label="pooled rank")],
              frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(f"{name.split('_')[0]} {label}: rank of the tissue's node per sample{ptxt}", x=0.5, y=0.97, fontsize=11.5)
    fig.savefig(out / "C_per_sample.png", dpi=180); fig.savefig(out / "C_per_sample.svg"); plt.close(fig)
    print(f"{name} ({tissue}) [linked{'' if a.score_key == 'scores' else ' ' + a.score_key}]: {verdict}{ptxt}  best target rank {best}/{L}: {summary['target_nodes'][:3]}  top: {summary['top_linked_nodes'][:4]}")


if __name__ == "__main__":
    main()

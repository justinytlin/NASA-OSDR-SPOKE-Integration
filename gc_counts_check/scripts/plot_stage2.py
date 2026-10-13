#!/usr/bin/env python3
"""Stage-2 summary figure: where each pilot tissue's own gene-linked node ranks, per SPOKE build.

Reads results/pilot_stage2_<tag>.tsv (pilot_stage2.py) for the tags given and the 2019 column they
carry, and draws one dot per build per tissue at rank / number of linked nodes (0 = top, as a
percentage). Tissues with no linked node on a build get an open marker at the right edge.
Usage: plot_stage2.py --tags 2025 2025noAeG 2025upAuG [2025nohk]
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; sys.path.insert(0, str(HERE)); import style as S

LABEL = {"2019": "2019 SPOKE (Zenodo PSEVs)", "2025": "2025 SPOKE, all Anatomy–Gene edges", "2025noAeG": "2025, without 'expresses' edges",
         "2025upAuG": "2025, 'upregulates' edges only", "2025nohk": "2025, all edges, housekeeping genes removed"}
COLOR = {"2019": "#2a78d6", "2025": "#eb6834", "2025noAeG": "#eda100", "2025upAuG": "#1baf7a", "2025nohk": "#4a3aa7"}
MARK = {"2019": "o", "2025": "s", "2025noAeG": "^", "2025upAuG": "D", "2025nohk": "v"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tags", nargs="+", default=["2025", "2025noAeG", "2025upAuG"])
    ap.add_argument("--no-2019", action="store_true", help="leave out the 2019 Zenodo series"); a = ap.parse_args()
    frames = {}
    for t in a.tags:
        d = pd.read_csv(ROOT / "results" / f"pilot_stage2_{t}.tsv", sep="\t")
        frames[t] = d.set_index("tissue")
    tissues = list(frames[a.tags[0]].index)
    base = frames[a.tags[0]]
    series = ([] if a.no_2019 else [("2019", base.rank_2019, base.n_linked_2019)]) + [(t, frames[t][f"rank_{t}"], frames[t][f"n_linked_{t}"]) for t in a.tags]
    S.setup()
    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(tissues) + 2.4)); h = 0.55 * len(tissues) + 2.4
    fig.subplots_adjust(left=0.2, right=0.97, top=1 - 0.6 / h, bottom=1.7 / h)
    yy = np.arange(len(tissues))[::-1]; off = np.linspace(-0.3, 0.3, len(series))
    for (tag, r, n), o in zip(series, off):
        pct = 100 * (r.values - 1) / (n.values - 1)
        ok = ~np.isnan(pct)
        ax.scatter(pct[ok], yy[ok] + o, s=46, color=COLOR[tag], marker=MARK[tag], edgecolors="white", linewidths=0.8, zorder=3)
        ax.scatter(np.full((~ok).sum(), 110), yy[~ok] + o, s=40, facecolors="white", edgecolors=COLOR[tag], marker=MARK[tag], linewidths=1.2, zorder=3)
    ax.axvline(0, color=S.AXIS, lw=0.8); ax.axvline(110, color=S.GRID, lw=0.8, ls=(0, (3, 3)))
    ax.set_yticks(yy); ax.set_yticklabels(tissues, fontsize=9.5); ax.set_xlim(-3, 114)
    ax.set_xticks([0, 25, 50, 75, 100, 110]); ax.set_xticklabels(["top", "25", "50", "75", "bottom", "no linked\nnode"], fontsize=9)
    ax.set_xlabel("rank of the tissue's own gene-linked Anatomy node among that build's linked nodes (% of the way down)", fontsize=9)
    ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], marker=MARK[t], ls="", color=COLOR[t], ms=7, label=LABEL[t]) for t, _, _ in series],
              loc="upper center", bbox_to_anchor=(0.45, -0.9 / (0.55 * len(tissues))), ncol=2, frameon=False, fontsize=8.5, handletextpad=0.5)
    ax.set_title("Own-tissue node rank by SPOKE 2025 edge set, ten pilot tissues (baseline-corrected)" if a.no_2019 else
                 "Own-tissue node rank by SPOKE build, ten pilot tissues (baseline-corrected)", fontsize=12, pad=10)
    out = ROOT / ("results/pilot_stage2_builds_2025only.png" if a.no_2019 else "results/pilot_stage2_builds.png"); fig.savefig(out, dpi=180); fig.savefig(out.with_suffix(".svg")); print("wrote", out)


if __name__ == "__main__":
    main()

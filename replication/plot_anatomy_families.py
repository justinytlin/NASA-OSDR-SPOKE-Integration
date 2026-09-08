#!/usr/bin/env python3
"""Figure: where do a tissue's Anatomy nodes sit in a PSEV rank table?

One row per rank table (e.g. the same dataset run on two PSEV builds). Left:
percentile-from-top of every Anatomy node, as strips for the named node
families (source tissue emphasised, controls in gray) over the distribution of
all Anatomy nodes. Right: the top-N Anatomy nodes with brain terms emphasised,
to show what actually occupies the extreme.

Usage:
  plot_anatomy_families.py --ranks FILE --label "2019 SPOKE" [--ranks FILE2 --label ...]
      [--col-regex RR10_HPC_GC]                columns to pool (default: all Rank_by_type_ columns)
      [--family "Hippocampal formation=hippocamp|dentate gyrus|Ammon|subiculum"]  (first = emphasised)
      [--family "Spleen=spleen|splenic"] [--family "Liver=liver|hepat"]
      [--brain-regex REGEX] [--top 12] [--out PATH.png]
"""
import argparse, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

# reference palette (dataviz skill): emphasis blue + text/chrome tokens
BLUE, INK, INK2, MUTED, GRID, AXIS, SURFACE = "#2a78d6", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
GRAY_MARK = "#b3b1aa"
STAGE_RE = re.compile(r"stage|life cycle|zygote|larval|pupal|neonate|embryo|fluid|secretion|milk|colostrum|"
                      r"processual|temporal boundary|organism|molecular system", re.I)


def load(path, col_regex):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cols = [c for c in df.columns if c.startswith("Rank_by_type_")]
    if col_regex:
        cols = [c for c in cols if re.search(col_regex, c, re.I)]
    an = df[df.Node_Type == "Anatomy"].copy()
    an["pct"] = 100 * an[cols].mean(axis=1) / len(an)
    return an, len(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", action="append", required=True)
    ap.add_argument("--label", action="append", default=[])
    ap.add_argument("--col-regex", default=None)
    ap.add_argument("--family", action="append", default=None)
    ap.add_argument("--brain-regex", default="brain|cortex|cerebr|hippocamp|telenceph|forebrain|gyrus|thalam|"
                                             "amygdala|striatum|pallium|dentate|nucleus|neur|spinal|ganglion|placode")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--title", default="Hippocampus in the ground-control expression PSEVs (OSD-564, 5 mice)")
    ap.add_argument("--out", default=str(Path(__file__).parent / "figures/anatomy_families.png"))
    a = ap.parse_args()
    fams = a.family or ["Hippocampal formation=hippocamp|dentate gyrus|Ammon|subiculum",
                        "Spleen=spleen|splenic", "Liver=liver|hepat"]
    fams = [f.split("=", 1) for f in fams]
    labels = a.label + [Path(p).stem for p in a.ranks[len(a.label):]]

    plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    rows = len(a.ranks)
    fig, axes = plt.subplots(rows, 2, figsize=(13, 3.6 * rows + 1.7), squeeze=False,
                             gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.75, "hspace": 0.6})
    fig.patch.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)

    for r, (path, lab) in enumerate(zip(a.ranks, labels)):
        an, ncols = load(path, a.col_regex)
        n = len(an)
        # ---------- left: family strips over the all-anatomy distribution ----------
        ax = axes[r, 0]; ax.set_facecolor(SURFACE)
        yrows = [f[0] for f in fams] + [f"All Anatomy ({n:,})"]
        ys = np.arange(len(yrows))[::-1]
        for i, (name, rx) in enumerate(fams):
            m = an.Node_Name.astype(str).str.contains(rx, case=False, regex=True)
            vals = an.pct[m].values
            emph = i == 0
            color = BLUE if emph else GRAY_MARK
            jitter = rng.uniform(-0.22, 0.22, len(vals))
            ax.scatter(vals, ys[i] + jitter, s=26 if emph else 18, color=color, alpha=0.9 if emph else 0.8,
                       edgecolors=SURFACE, linewidths=1.2, zorder=3)
            med = np.median(vals)
            ax.plot([med, med], [ys[i] - 0.36, ys[i] + 0.36], color=INK if emph else INK2, lw=2.2, zorder=4)
            p = mannwhitneyu(vals, an.pct[~m], alternative="less").pvalue
            ptxt = "p < 1e-15" if p < 1e-15 else f"p = {p:.2g}"
            yrows[i] = f"{name} ({m.sum()} nodes)\nmedian {med:.0f}%, {ptxt}"
        # all-anatomy row: histogram strip
        h, edges = np.histogram(an.pct, bins=50, range=(0, 100))
        ax.bar((edges[:-1] + edges[1:]) / 2, h / h.max() * 0.7, width=2.0 * 0.92, bottom=ys[-1] - 0.35,
               color=GRID, edgecolor="none", zorder=2)
        ax.plot([50, 50], [ys[-1] - 0.36, ys[-1] + 0.36], color=INK2, lw=2.2, zorder=4)
        yrows[-1] = f"All Anatomy ({n:,} nodes)\nmedian 50%"
        # top 2.5% marker
        ax.axvline(2.5, color=INK2, lw=0.8, zorder=1)
        ax.text(2.5, len(yrows) - 0.42, "top 2.5%", ha="left", va="bottom", fontsize=8.5, color=INK2)
        ax.set_yticks(ys); ax.set_yticklabels(yrows, fontsize=9.5, linespacing=1.4)
        ax.set_xlim(-1, 101); ax.set_ylim(-0.6, len(yrows) - 0.35)
        ax.set_xlabel("Anatomy rank percentile within the sample profiles (0 = most associated with expressed genes)",
                      color=INK2, fontsize=9)
        ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xticklabels(["0", "25", "50", "75", "100%"])
        ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
        for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_title(f"{lab}: {ncols} samples pooled", loc="left", fontsize=11, color=INK, pad=12)

        # ---------- right: what occupies the top ----------
        ax2 = axes[r, 1]; ax2.set_facecolor(SURFACE)
        top = an.nsmallest(a.top, "pct")
        names = top.Node_Name.astype(str).tolist()
        is_brain = [bool(re.search(a.brain_regex, s, re.I)) and not STAGE_RE.search(s) for s in names]
        yy = np.arange(len(names))[::-1]
        ranks = top.pct.values / 100 * n
        ax2.hlines(yy, 0, ranks, color=GRID, lw=1.2, zorder=1)
        ax2.scatter(ranks, yy, s=[46 if b else 30 for b in is_brain],
                    color=[BLUE if b else GRAY_MARK for b in is_brain], edgecolors=SURFACE, linewidths=1.2, zorder=3)
        ax2.set_yticks(yy)
        ax2.set_yticklabels([("* " if STAGE_RE.search(s) else "") + s for s in names], fontsize=8.5)
        ax2.set_xlim(0, ranks.max() * 1.15)
        ax2.set_xlabel(f"mean rank among {n:,} Anatomy nodes (1 = top)", color=INK2, fontsize=9)
        ax2.grid(axis="x", color=GRID, lw=0.8); ax2.set_axisbelow(True)
        for s in ("top", "right", "left"): ax2.spines[s].set_visible(False)
        ax2.tick_params(axis="y", length=0)
        nb = sum(is_brain)
        ax2.set_title(f"What sits at the top: {nb} of {a.top} are brain/neural", loc="left",
                      fontsize=11, color=INK, pad=12)

    handles = [Line2D([], [], marker="o", ls="", color=BLUE, markersize=7, label=f"{fams[0][0]} / brain & neural terms"),
               Line2D([], [], marker="o", ls="", color=GRAY_MARK, markersize=6, label="control tissues / other anatomy"),
               Line2D([], [], color=INK, lw=2.2, label="median")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.03, 0.01), ncol=3, frameon=False,
               fontsize=9, handletextpad=0.5, columnspacing=1.6)
    fig.suptitle(a.title, x=0.03, ha="left", fontsize=13.5, color=INK, y=0.985, fontweight="semibold")
    fig.text(0.03, 0.925, "Each dot is one Anatomy node, placed by its mean rank-by-type across the samples as a "
             "percentile of all Anatomy nodes.   * = UBERON life-stage / organism pseudo-anatomy",
             fontsize=9, color=INK2)
    top_frac = 1 - 1.15 / (3.6 * rows + 1.7)
    bot_frac = 0.85 / (3.6 * rows + 1.7)
    fig.subplots_adjust(top=top_frac, bottom=bot_frac, left=0.2, right=0.985)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, facecolor=SURFACE)
    fig.savefig(out.with_suffix(".svg"), facecolor=SURFACE)
    print(f"wrote {out} and {out.with_suffix('.svg')}")


if __name__ == "__main__":
    main()

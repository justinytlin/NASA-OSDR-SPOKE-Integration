"""Figure for the cell-type route: rank of the tissue's matching CellType node for datasets whose
tissue has no gene-linked Anatomy node in a build.

usage: python scripts/plot_celltype.py --tag 2025upAuG [--all]
reads  results/celltype_check_<tag>.tsv
writes results/celltype_nolinked_<tag>.png/.svg  (or celltype_all_<tag> with --all)
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import style as S

ROOT = Path(__file__).resolve().parents[1]
NO_LINKED = ["retina", "bone marrow", "placenta", "cecum", "thymus", "eye", "optic nerve", "bone", "mandible", "temporal bone"]
VERDICT = {"CLEAR": ("#0ca30c", "o"), "PARTIAL": ("#fab219", "^"), "NONE": ("#d03b3b", "X")}


def verdict(r):
    return "CLEAR" if r == 1 else "PARTIAL" if r <= 7 else "NONE"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    df = pd.read_csv(ROOT / "results" / f"celltype_check_{a.tag}.tsv", sep="\t")
    df = df[df.celltype_rank.notna()].copy(); df["celltype_rank"] = df.celltype_rank.astype(int)
    if not a.all:
        df = df[df.tissue.isin(NO_LINKED)]
    df["verdict"] = df.celltype_rank.map(verdict)
    # order tissues by median rank, datasets within tissue by rank
    med = df.groupby("tissue").celltype_rank.median().sort_values()
    df["tkey"] = df.tissue.map({t: i for i, t in enumerate(med.index)})
    df = df.sort_values(["tkey", "celltype_rank", "dataset"]).reset_index(drop=True)
    n_ct = int(df.n_celltypes.iloc[0])
    S.setup(); n = len(df)
    # layout rows: a header row per tissue, then one row per dataset
    rows = []  # (kind, payload)
    for t in med.index:
        rows.append(("head", t))
        for _, r in df[df.tissue == t].iterrows():
            rows.append(("data", r))
    nr = len(rows); h = 0.34 * nr + 2.4
    fig, ax = plt.subplots(figsize=(10, h)); fig.subplots_adjust(left=0.22, right=0.98, top=1 - 0.6 / h, bottom=1.35 / h)
    yy = np.arange(nr)[::-1]
    ticks_y, ticks_l = [], []
    for y, (kind, r) in zip(yy, rows):
        if kind == "head":
            ax.text(0.86, y, r, ha="right", va="center", fontsize=10, color=S.INK, fontweight="bold", clip_on=False)
            if y != yy[0]:
                ax.axhline(y + 0.5, color=S.GRID, lw=0.8)
            continue
        c, mk = VERDICT[r.verdict]
        ax.scatter(r.celltype_rank, y, s=52, marker=mk, color=c, edgecolors="white", linewidths=0.8, zorder=3)
        ax.text(r.celltype_rank * 1.25, y, r.best_family_celltype, va="center", ha="left", fontsize=8.5, color=S.INK2)
        ticks_y.append(y); ticks_l.append(f"{r.dataset.replace('_', ' ')}  ({r.n_samples})")
    ax.set_yticks(ticks_y); ax.set_yticklabels(ticks_l, fontsize=9)
    ax.set_xscale("log"); ax.set_xlim(0.8, n_ct * 1.15)
    ticks = [1, 2, 5, 10, 20, 50, 100, 500, n_ct]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
    ax.axvline(7.5, color=S.GRID, lw=0.8, ls=(0, (3, 3)))
    ax.set_xlabel(f"rank of the tissue's matching CellType node among {n_ct:,} cell types (1 = top)")
    ax.set_ylim(-0.7, nr - 0.3)
    ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
    scope = "all tissues" if a.all else "tissues without a linked Anatomy node"
    ax.set_title(f"Cell-type route, {a.tag}: {scope} ({n} datasets)", loc="left", fontsize=12, pad=12)
    ax.legend(handles=[Line2D([], [], marker=mk, ls="", color=c, ms=7, label=f"{v.lower()} ({'rank 1' if v=='CLEAR' else 'rank 2–7' if v=='PARTIAL' else 'below 7'})")
                       for v, (c, mk) in VERDICT.items()],
              loc="upper center", bbox_to_anchor=(0.5, -0.75 / h - 0.03), ncol=3, frameon=False, fontsize=9)
    stem = ROOT / "results" / (f"celltype_all_{a.tag}" if a.all else f"celltype_nolinked_{a.tag}")
    fig.savefig(f"{stem}.png", dpi=170); fig.savefig(f"{stem}.svg"); print("wrote", stem)


if __name__ == "__main__":
    main()

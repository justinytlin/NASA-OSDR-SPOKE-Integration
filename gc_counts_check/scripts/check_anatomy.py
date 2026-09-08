#!/usr/bin/env python3
"""Single-dataset check: does the source tissue separate in the SPOKE Anatomy ranking?

Uses ONLY the given dataset (no comparison with other datasets). Input is a
per-sample rank table from gc_counts_to_psev.py / replication/counts_to_psev.py
(any table with Rank_by_type_* columns works, e.g. a Log2fc DE-route table), or
a GeneLab unnormalized counts CSV, in which case the 2019-PSEV projection is run
first (~7 min).

Evidence (all within the dataset):
  QQ            observed percentiles of the tissue family's Anatomy nodes vs the
                percentiles a random set of Anatomy nodes would have (diagonal);
                95 % envelope from 2,000 random node sets of the same size;
                one-sided KS test and Mann-Whitney U vs all other Anatomy nodes
  competition   the tissue family's median percentile against the same statistic
                for every other organ family in this dataset (rank among families)
  consistency   per-sample median percentile and per-sample Mann-Whitney p

Verdict (printed and written to <out>/<name>/summary.json):
  CLEAR    pooled MWU p < 0.01, KS p < 0.01, >= 80 % of samples individually
           p < 0.05, and the family is the best-ranked of all organ families
  PARTIAL  pooled MWU p < 0.05 but one or more of the other criteria fail
  NONE     pooled MWU p >= 0.05
Thresholds are flags (--p-pooled, --p-sample, --min-sample-frac).

Usage:
  check_anatomy.py --ranks TABLE.tsv --tissue spleen [--name OSD-246]
  check_anatomy.py --counts GLDS-xxx_..._Unnormalized_Counts_GLbulkRNAseq.csv --tissue liver --name OSD-245 [--samples _GC_]
  --tissue accepts a built-in family name (hippocampus, thymus, spleen, liver, heart, kidney,
  lung, skeletal muscle, bone marrow, skin, intestine, testis, blood, cerebral cortex, adrenal gland,
  eye, retina, brain) or
  --family-regex / --exclude-regex for anything else.
Output folder <out>/<name>/ (default out = GC Counts Check/datasets):
  A_qq_plot.png/.svg, B_organ_families.png/.svg, C_per_sample.png/.svg  one figure each
  summary.json (verdict + criteria + numbers), organ_families.tsv, per_sample.tsv
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import style as S
from families import FAMILIES, CORE_ORGAN_FAMILIES

STAGE_RE = re.compile(r"stage|life cycle|zygote|larval|pupal|neonate|embryo|fluid|secretion|milk|colostrum|"
                      r"processual|temporal boundary|organism|molecular system", re.I)
ORGAN_FAMILIES = list(CORE_ORGAN_FAMILIES)


def family_mask(names, inc, exc):
    m = names.str.contains(inc, case=False, regex=True)
    if exc:
        m &= ~names.str.contains(exc, case=False, regex=True)
    return m.values


def run_projection(counts, name, tissue, samples, threads):
    cmd = [sys.executable, str(HERE / "gc_counts_to_psev.py"), "--counts", str(counts), "--name", name,
           "--tissue", tissue, "--threads", str(threads)]
    for s in samples or []:
        cmd += ["--samples", s]
    print("running projection:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return ROOT / "psev_out" / f"{name}_GC_counts_ranks_and_rank_by_type_for_meta0_1.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranks", help="rank table with Rank_by_type_* columns")
    ap.add_argument("--counts", help="GeneLab unnormalized counts CSV (runs the projection first)")
    ap.add_argument("--tissue", required=True, help="source tissue (built-in family name or free label with --family-regex)")
    ap.add_argument("--family-regex", default=None); ap.add_argument("--exclude-regex", default=None)
    ap.add_argument("--name", default=None); ap.add_argument("--samples", action="append", default=None)
    ap.add_argument("--col-regex", default=None, help="restrict to Rank_by_type_ columns matching this")
    ap.add_argument("--out", default=str(ROOT / "datasets"))
    ap.add_argument("--p-pooled", type=float, default=0.01); ap.add_argument("--p-sample", type=float, default=0.05)
    ap.add_argument("--min-sample-frac", type=float, default=0.8)
    ap.add_argument("--n-null", type=int, default=2000); ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    if not a.ranks and not a.counts:
        sys.exit("give --ranks or --counts")
    name = a.name or (Path(a.ranks or a.counts).stem.split("_")[0])
    ranks_path = Path(a.ranks) if a.ranks else run_projection(a.counts, name, a.tissue, a.samples, a.threads)
    tissue = a.tissue
    if a.family_regex:
        inc, exc = a.family_regex, a.exclude_regex
    elif tissue in FAMILIES:
        inc, exc = FAMILIES[tissue]
    else:
        sys.exit(f"'{tissue}' is not a built-in family ({', '.join(ORGAN_FAMILIES)}); pass --family-regex")
    color = S.TISSUE_COLOR.get(tissue, "#2a78d6")
    out = Path(a.out) / name; out.mkdir(parents=True, exist_ok=True)

    # ---------- load ----------
    df = pd.read_csv(ranks_path, sep="\t", low_memory=False)
    cols = [c for c in df.columns if c.startswith("Rank_by_type_")]
    if a.col_regex:
        cols = [c for c in cols if re.search(a.col_regex, c, re.I)]
    an = df[df.Node_Type == "Anatomy"].reset_index(drop=True)
    names = an.Node_Name.astype(str)
    n = len(an)
    pct = an[cols].values * 100.0 / n                    # nodes x samples
    pooled = pct.mean(axis=1)
    samples = [c[len("Rank_by_type_"):] for c in cols]
    m = family_mask(names, inc, exc)
    if m.sum() < 5:
        sys.exit(f"only {m.sum()} Anatomy nodes match the family; widen the regex")
    rng = np.random.default_rng(0)

    # ---------- statistics ----------
    fam, rest = pooled[m], pooled[~m]
    mwu = mannwhitneyu(fam, rest, alternative="less").pvalue
    ks = ks_2samp(fam, rest, alternative="greater").pvalue
    per = pd.DataFrame({"sample": samples,
                        "family_median_pct": [np.median(pct[m, j]) for j in range(len(samples))],
                        "mwu_p_less": [mannwhitneyu(pct[m, j], pct[~m, j], alternative="less").pvalue for j in range(len(samples))]})
    frac_sig = float((per.mwu_p_less < a.p_sample).mean())
    comp = []
    for f in ORGAN_FAMILIES:
        mf = family_mask(names, *FAMILIES[f])
        if mf.sum() >= 5:
            comp.append(dict(family=f, n_nodes=int(mf.sum()), median_pct=float(np.median(pooled[mf])),
                             mwu_p_less=float(mannwhitneyu(pooled[mf], pooled[~mf], alternative="less").pvalue)))
    if tissue not in [c["family"] for c in comp]:
        comp.append(dict(family=tissue, n_nodes=int(m.sum()), median_pct=float(np.median(fam)), mwu_p_less=float(mwu)))
    comp = pd.DataFrame(comp).sort_values("median_pct").reset_index(drop=True)
    fam_rank = int(comp.index[comp.family == tissue][0]) + 1
    top_real = an.loc[~names.str.contains(STAGE_RE)].assign(pct=pooled[~names.str.contains(STAGE_RE).values]).nsmallest(10, "pct")
    top_fam = an.loc[m].assign(pct=pooled[m]).nsmallest(5, "pct")

    if mwu < a.p_pooled and ks < a.p_pooled and frac_sig >= a.min_sample_frac and fam_rank == 1:
        verdict = "CLEAR"
    elif mwu < 0.05:
        verdict = "PARTIAL"
    else:
        verdict = "NONE"
    reasons = [f"pooled Mann-Whitney p = {mwu:.2g} ({'ok' if mwu < a.p_pooled else 'fail'})",
               f"KS p = {ks:.2g} ({'ok' if ks < a.p_pooled else 'fail'})",
               f"{int((per.mwu_p_less < a.p_sample).sum())}/{len(per)} samples individually p < {a.p_sample} "
               f"({'ok' if frac_sig >= a.min_sample_frac else 'fail'})",
               f"family ranks {fam_rank} of {len(comp)} organ families by median percentile ({'ok' if fam_rank == 1 else 'fail'})"]
    summary = dict(name=name, tissue=tissue, ranks=str(ranks_path), n_samples=len(samples), n_anatomy=n,
                   family_regex=inc, family_exclude=exc, n_family_nodes=int(m.sum()), family_median_pct=float(np.median(fam)),
                   in_top_2_5pct=int((fam <= 2.5).sum()), in_top_10pct=int((fam <= 10).sum()),
                   mwu_p_less=float(mwu), ks_p_greater=float(ks), frac_samples_significant=frac_sig,
                   family_rank_among_organ_families=fam_rank, verdict=verdict, criteria=reasons,
                   best_family_nodes=[f"{r.Node_Name} ({r.pct:.1f}%)" for r in top_fam.itertuples()],
                   top_real_anatomy=[f"{r.Node_Name} ({r.pct:.2f}%)" for r in top_real.itertuples()])
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    comp.to_csv(out / "organ_families.tsv", sep="\t", index=False)
    per.to_csv(out / "per_sample.tsv", sep="\t", index=False)
    verdict_line = (f"Verdict {verdict}: pooled Mann-Whitney p = {mwu:.2g}, KS p = {ks:.2g}, "
                    f"{int((per.mwu_p_less < a.p_sample).sum())}/{len(per)} samples p < {a.p_sample}, "
                    f"family ranks {fam_rank} of {len(comp)} organ families")

    def finish(fig, letter, title, sub, fname):
        fig.suptitle(f"{name} {tissue}: {title}", x=0.5, ha="center", y=0.97, fontsize=12.5)
        fig.savefig(out / f"{fname}.png", dpi=180); fig.savefig(out / f"{fname}.svg"); plt.close(fig)

    # ---------- figures (one file each) ----------
    S.setup()
    k = int(m.sum()); exp = 100.0 * (np.arange(1, k + 1) - 0.5) / k
    n_s = len(samples)

    # A. QQ vs uniform
    fig, ax = plt.subplots(figsize=(6.4, 6.6)); fig.subplots_adjust(left=0.13, right=0.97, top=0.92, bottom=0.2)
    draws = np.sort(np.array([rng.choice(pooled, k, replace=False) for _ in range(a.n_null)]), axis=1)
    ax.fill_between(exp, np.percentile(draws, 2.5, axis=0), np.percentile(draws, 97.5, axis=0), color=S.GRID, alpha=0.85, lw=0, zorder=1)
    ax.plot([0, 100], [0, 100], color=S.AXIS, lw=1, ls=(0, (4, 3)), zorder=2)
    for j in range(n_s):
        ax.plot(exp, np.sort(pct[m, j]), color=color, lw=0.8, alpha=0.3, zorder=3)
    for f in ORGAN_FAMILIES:
        if f == tissue:
            continue
        mf = family_mask(names, *FAMILIES[f]); kf = int(mf.sum())
        if kf >= 5:
            ax.plot(100.0 * (np.arange(1, kf + 1) - 0.5) / kf, np.sort(pooled[mf]), color=S.GRAY_MARK, lw=0.9, alpha=0.7, zorder=3)
    ax.plot(exp, np.sort(fam), color=color, lw=2.4, marker="o", ms=3.5, mec=S.SURFACE, mew=0.6, zorder=5)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_aspect("equal")
    ax.set_xlabel("expected percentile if the family were random Anatomy nodes", fontsize=9)
    ax.set_ylabel("observed percentile of the family's nodes (0 = top)", fontsize=9)
    ax.grid(color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax, keep=("bottom", "left"))
    handles = [Line2D([], [], color=color, lw=2.4, marker="o", ms=4, label=f"{tissue} family, samples pooled"),
               Line2D([], [], color=color, lw=0.8, alpha=0.5, label="each GC sample"),
               Line2D([], [], color=S.GRAY_MARK, lw=0.9, label="other organ families"),
               plt.Rectangle((0, 0), 1, 1, color=S.GRID, label="95 % envelope, random node sets"),
               Line2D([], [], color=S.AXIS, lw=1, ls=(0, (4, 3)), label="no association")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.005), ncol=2, frameon=False, fontsize=8.5)
    finish(fig, "A", f"QQ plot of the {tissue} Anatomy family ({k} nodes)", "", "A_qq_plot")

    # B. competition among organ families
    fig, ax = plt.subplots(figsize=(8.2, 6.0)); fig.subplots_adjust(left=0.22, right=0.97, top=0.92, bottom=0.11)
    yy = np.arange(len(comp))[::-1]
    for i, r in comp.iterrows():
        mf = family_mask(names, *FAMILIES[r.family]) if r.family in FAMILIES else m
        vals = pooled[mf]; own = r.family == tissue
        ax.scatter(vals, yy[i] + rng.uniform(-0.2, 0.2, len(vals)), s=16 if own else 9, color=color if own else S.GRAY_MARK,
                   alpha=0.85 if own else 0.55, edgecolors=S.SURFACE, linewidths=0.6, zorder=3)
        ax.plot([r.median_pct] * 2, [yy[i] - 0.38, yy[i] + 0.38], color=S.INK if own else S.INK2, lw=2.2 if own else 1.6, zorder=4)
    ax.axvline(50, color=S.AXIS, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels([f"{r.family} ({r.n_nodes})" + ("  ◀" if r.family == tissue else "") for _, r in comp.iterrows()],
                                          fontsize=8.5)
    for lab, (_, r) in zip(ax.get_yticklabels(), comp.iterrows()):
        if r.family == tissue:
            lab.set_color(S.INK); lab.set_fontweight("semibold")
    ax.set_xlim(0, 100); ax.set_xlabel("pooled percentile of each family's Anatomy nodes (black bar = family median; 0 = top)", fontsize=9)
    ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
    finish(fig, "B", f"organ families ({tissue} ranks {fam_rank} of {len(comp)})", "", "B_organ_families")

    # C. per-sample consistency
    fig, ax = plt.subplots(figsize=(max(6.4, 0.32 * n_s + 3.2), 5.6)); fig.subplots_adjust(left=0.16, right=0.97, top=0.9, bottom=0.3)
    xs = np.arange(n_s); sig = per.mwu_p_less.values < a.p_sample
    ax.axhline(50, color=S.AXIS, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.axhline(np.median(fam), color=S.INK2, lw=1, zorder=2)
    ax.scatter(xs, per.family_median_pct, s=[46 if s_ else 30 for s_ in sig], color=[color if s_ else S.SURFACE for s_ in sig],
               edgecolors=color, linewidths=1.4, zorder=3)
    ax.set_xticks(xs); ax.set_xticklabels([re.sub(r"^(Mmus_C57-6[A-Z]*_[A-Z]+_|RR10_[A-Z]+_)", "", s_)[:22] for s_ in samples], rotation=90, fontsize=7)
    ax.set_ylim(0, 100); ax.invert_yaxis()
    ax.set_ylabel("family median percentile in the sample (0 = top)", fontsize=9)
    ax.grid(axis="y", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax, keep=("bottom", "left"))
    handles = [Line2D([], [], marker="o", ls="", color=color, mec=color, ms=7, label=f"sample p < {a.p_sample}"),
               Line2D([], [], marker="o", ls="", color=S.SURFACE, mec=color, mew=1.4, ms=6, label="not significant"),
               Line2D([], [], color=S.INK2, lw=1, label="pooled median")]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    finish(fig, "C", f"{tissue} family per sample ({int(sig.sum())}/{n_s} samples p < {a.p_sample})", "", "C_per_sample")

    print(f"\n{name} ({tissue}): {verdict}")
    for r in reasons:
        print("  -", r)
    print(f"  best family nodes: {', '.join(summary['best_family_nodes'])}")
    print(f"  top real anatomy nodes (hubs, for context): {', '.join(summary['top_real_anatomy'][:5])}")
    print(f"  wrote {out}/ (A_qq_plot, B_organ_families, C_per_sample .png/.svg; summary.json; organ_families.tsv; per_sample.tsv)")


if __name__ == "__main__":
    main()

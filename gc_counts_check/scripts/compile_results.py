#!/usr/bin/env python3
"""Compile every dataset's anatomy-check summary into one table and a per-tissue overview.

Each dataset was judged on its own (datasets/<name>/summary.json); this only tabulates
those verdicts. Outputs in results/:
  compiled_verdicts.tsv          one row per dataset (study x tissue)
  tissue_summary.tsv             one row per tissue: number of datasets and CLEAR / PARTIAL / NONE counts
  tissue_separation_overview.png each dataset's family median percentile, grouped by tissue, marked by verdict
"""
import json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import style as S

LEGACY = {  # first two sets, run before the manifest existed
    "OSD-564": ("RR-10", "hippocampus"), "OSD-244": ("RR-6", "thymus"), "OSD-245": ("RR-6", "liver"), "OSD-246": ("RR-6", "spleen"),
    "OSD-288": ("MHU-1", "spleen"), "OSD-289": ("MHU-1, MHU-2", "thymus"), "OSD-137": ("RR-3", "liver"),
    "OSD-161": ("RR-3", "adrenal gland"), "OSD-162": ("RR-3", "eye"), "OSD-163": ("RR-3", "kidney"), "OSD-168": ("RR-3", "liver"),
    "OSD-194": ("RR-3", "retina"), "OSD-270": ("RR-3", "heart"), "OSD-352": ("RR-3", "brain")}
VERDICT = {"CLEAR": ("#0ca30c", "o"), "PARTIAL": ("#fab219", "^"), "NONE": ("#d03b3b", "X")}


def main():
    mf = ROOT / "data/_catalog/manifest.json"
    man = {e["name"]: e for e in json.load(open(mf))} if mf.exists() else {}
    rows = []
    for sj in sorted((ROOT / "datasets").glob("*/summary.json")):
        s = json.load(open(sj)); name = s["name"]; e = man.get(name)
        fam_tab = pd.read_csv(sj.parent / "organ_families.tsv", sep="\t")
        payload, label = (e["payload"], e["tissue"]) if e else LEGACY.get(name, ("", s["tissue"]))
        rows.append(dict(dataset=name, accession=name.split("_")[0], payload=payload, tissue=label, family=s["tissue"],
                         control_type=e["control_type"] if e else "ground control", n_samples=s["n_samples"],
                         n_family_nodes=s["n_family_nodes"], family_median_pct=round(s["family_median_pct"], 2),
                         in_top_2_5pct=s["in_top_2_5pct"], in_top_10pct=s["in_top_10pct"], mwu_p=s["mwu_p_less"], ks_p=s["ks_p_greater"],
                         frac_samples_significant=round(s["frac_samples_significant"], 2),
                         family_rank=s["family_rank_among_organ_families"], n_organ_families=len(fam_tab), verdict=s["verdict"],
                         best_family_nodes="; ".join(s["best_family_nodes"][:3])))
        cj = sj.parent / "baseline_corrected" / "summary.json"
        if cj.exists():
            c = json.load(open(cj)); ctab = pd.read_csv(cj.parent / "organ_families.tsv", sep="\t")
            rows[-1].update(corr_family_median_pct=round(c["family_median_pct"], 2), corr_in_top_2_5pct=c["in_top_2_5pct"],
                            corr_in_top_10pct=c["in_top_10pct"], corr_mwu_p=c["mwu_p_less"], corr_ks_p=c["ks_p_greater"],
                            corr_frac_samples_significant=round(c["frac_samples_significant"], 2),
                            corr_family_rank=c["family_rank_among_organ_families"], corr_n_organ_families=len(ctab),
                            corr_verdict=c["verdict"], corr_best_family_nodes="; ".join(c["best_family_nodes"][:3]),
                            corr_top_anatomy="; ".join(c["top_real_anatomy"][:5]))
        oj = sj.parent / "corrected_ontology" / "summary.json"
        if oj.exists():
            c = json.load(open(oj)); otab = pd.read_csv(oj.parent / "organ_families.tsv", sep="\t")
            rows[-1].update(onto_n_family_nodes=c["n_family_nodes"], onto_family_median_pct=round(c["family_median_pct"], 2),
                            onto_mwu_p=c["mwu_p_less"], onto_ks_p=c["ks_p_greater"],
                            onto_frac_samples_significant=round(c["frac_samples_significant"], 2),
                            onto_family_rank=c["family_rank_among_organ_families"], onto_n_organ_families=len(otab),
                            onto_verdict=c["verdict"], onto_best_family_nodes="; ".join(c["best_family_nodes"][:3]))
        lj = sj.parent / "corrected_linked" / "summary.json"
        if lj.exists():
            c = json.load(open(lj))
            rows[-1].update(linked_verdict=c["verdict"] + (" (proxy)" if c.get("proxy") and c["verdict"] != "NO LINKED NODE" else ""),
                            linked_best_target_rank=c.get("best_target_rank"), linked_n=c.get("n_linked"),
                            linked_targets="; ".join(c.get("target_nodes", [])[:3]),
                            linked_frac_samples_rank1=c.get("frac_samples_target_rank1"),
                            linked_top_nodes="; ".join(c.get("top_linked_nodes", [])[:5]))
    df = pd.DataFrame(rows)
    out = ROOT / "results"; out.mkdir(exist_ok=True)
    order = {"CLEAR": 0, "PARTIAL": 1, "NONE": 2}
    df = df.sort_values(["tissue", "family_median_pct"]).reset_index(drop=True)
    df.to_csv(out / "compiled_verdicts.tsv", sep="\t", index=False)
    ts = df.groupby("tissue").agg(family=("family", "first"), n_datasets=("dataset", "size"), n_samples=("n_samples", "sum"),
                                  CLEAR=("verdict", lambda v: int((v == "CLEAR").sum())), PARTIAL=("verdict", lambda v: int((v == "PARTIAL").sum())),
                                  NONE=("verdict", lambda v: int((v == "NONE").sum())),
                                  median_family_pct=("family_median_pct", "median"), min_family_pct=("family_median_pct", "min"),
                                  max_family_pct=("family_median_pct", "max"),
                                  datasets_pooled_p_lt_0_01=("mwu_p", lambda p: int((p < 0.01).sum()))).reset_index()
    ts["separates"] = np.where(ts.CLEAR + ts.PARTIAL == ts.n_datasets, "yes", np.where(ts.CLEAR + ts.PARTIAL == 0, "no", "mixed"))
    has_corr = "corr_verdict" in df and df.corr_verdict.notna().any()
    if has_corr:
        dc = df[df.corr_verdict.notna()]
        tc = dc.groupby("tissue").agg(corr_n=("dataset", "size"), corr_CLEAR=("corr_verdict", lambda v: int((v == "CLEAR").sum())),
                                      corr_PARTIAL=("corr_verdict", lambda v: int((v == "PARTIAL").sum())),
                                      corr_NONE=("corr_verdict", lambda v: int((v == "NONE").sum())),
                                      corr_median_family_pct=("corr_family_median_pct", "median"),
                                      corr_median_family_rank=("corr_family_rank", "median")).reset_index()
        tc["corr_separates"] = np.where(tc.corr_CLEAR + tc.corr_PARTIAL == tc.corr_n, "yes",
                                        np.where(tc.corr_CLEAR + tc.corr_PARTIAL == 0, "no", "mixed"))
        ts = ts.merge(tc, on="tissue", how="left")
    if "onto_verdict" in df and df.onto_verdict.notna().any():
        to = df[df.onto_verdict.notna()].groupby("tissue").agg(
            onto_CLEAR=("onto_verdict", lambda v: int((v == "CLEAR").sum())), onto_PARTIAL=("onto_verdict", lambda v: int((v == "PARTIAL").sum())),
            onto_NONE=("onto_verdict", lambda v: int((v == "NONE").sum())), onto_family_nodes=("onto_n_family_nodes", "first"),
            onto_median_family_pct=("onto_family_median_pct", "median"), onto_median_family_rank=("onto_family_rank", "median")).reset_index()
        ts = ts.merge(to, on="tissue", how="left")
    if "linked_verdict" in df and df.linked_verdict.notna().any():
        tl = df[df.linked_verdict.notna()].groupby("tissue").agg(
            linked_target=("linked_targets", lambda v: re.sub(r" \(rank \d+\)", "", str(v.iloc[0])) if isinstance(v.iloc[0], str) else ""),
            linked_proxy=("linked_verdict", lambda v: bool(v.str.contains("proxy").any())),
            linked_median_rank=("linked_best_target_rank", "median"),
            linked_CLEAR=("linked_verdict", lambda v: int(v.str.startswith("CLEAR").sum())),
            linked_PARTIAL=("linked_verdict", lambda v: int(v.str.startswith("PARTIAL").sum())),
            linked_NONE=("linked_verdict", lambda v: int(v.str.startswith("NONE").sum())),
            linked_no_node=("linked_verdict", lambda v: int(v.str.startswith("NO LINKED").sum()))).reset_index()
        ts = ts.merge(tl, on="tissue", how="left")
    ts = ts.sort_values("corr_median_family_pct" if has_corr else "median_family_pct").reset_index(drop=True)
    ts.to_csv(out / "tissue_summary.tsv", sep="\t", index=False)

    S.setup()

    def overview(pcol, vcol, title, fname):
        d0 = df[df[vcol].notna()]
        order = d0.groupby("tissue")[pcol].median().sort_values().index.tolist()
        h = 0.33 * len(order) + 1.8
        fig, ax = plt.subplots(figsize=(9, h)); fig.subplots_adjust(left=0.27, right=0.97, top=1 - 0.5 / h, bottom=1.0 / h)
        yy = np.arange(len(order))[::-1]; rng = np.random.default_rng(0)
        for y, t in zip(yy, order):
            d = d0[d0.tissue == t]
            for v, (c, mk) in VERDICT.items():
                g = d[d[vcol] == v]
                ax.scatter(g[pcol], y + rng.uniform(-0.16, 0.16, len(g)), s=46, color=c, marker=mk, edgecolors="white", linewidths=0.8, zorder=3)
        ax.axvline(50, color=S.AXIS, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.set_yticks(yy); ax.set_yticklabels([f"{t} ({(d0.tissue == t).sum()})" for t in order], fontsize=9)
        ax.set_xlim(0, 100); ax.set_ylim(-0.7, len(order) - 0.3)
        ax.set_xlabel("median percentile of the tissue's Anatomy family in the dataset (0 = top, 50 = no association)", fontsize=9)
        ax.grid(axis="x", color=S.GRID, lw=0.7); ax.set_axisbelow(True); S.despine(ax); ax.tick_params(axis="y", length=0)
        ax.legend(handles=[Line2D([], [], marker=mk, ls="", color=c, ms=7, label=v) for v, (c, mk) in VERDICT.items()],
                  loc="upper center", bbox_to_anchor=(0.5, -0.9 / (0.33 * len(order))), ncol=3, frameon=False, fontsize=9)
        ax.set_title(f"{title}: {len(d0)} datasets, {len(order)} tissues", fontsize=12.5, pad=10)
        fig.savefig(out / f"{fname}.png", dpi=180); fig.savefig(out / f"{fname}.svg"); plt.close(fig)

    overview("family_median_pct", "verdict", "Tissue separation in control samples (raw score)", "tissue_separation_overview")
    if has_corr:
        overview("corr_family_median_pct", "corr_verdict", "Tissue separation in control samples (baseline-corrected)",
                 "tissue_separation_overview_corrected")
    if "onto_verdict" in df and df.onto_verdict.notna().any():
        overview("onto_family_median_pct", "onto_verdict", "Tissue separation (baseline-corrected, ontology families)",
                 "tissue_separation_overview_ontology")
    print(f"{len(df)} datasets, {len(ts)} tissues;", df.verdict.value_counts().to_dict())
    cols = ["tissue", "n_datasets", "CLEAR", "PARTIAL", "NONE", "median_family_pct", "separates"]
    if has_corr:
        cols += ["corr_n", "corr_CLEAR", "corr_PARTIAL", "corr_NONE", "corr_median_family_pct", "corr_separates"]
        print("corrected:", df.corr_verdict.value_counts().to_dict())
    print(ts[cols].to_string(index=False))


if __name__ == "__main__":
    main()

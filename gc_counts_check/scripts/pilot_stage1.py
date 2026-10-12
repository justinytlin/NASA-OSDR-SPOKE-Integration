#!/usr/bin/env python3
"""Stage-1 pilot: four ways of scoring, ten tissues, judged on the gene-linked nodes and on ontology families.

Arms (all baseline-corrected, each dataset on its own):
  A  rank-based projection (the current method)                 <base>__pilot  score-key scores
  Z  z-score projection (no within-gene rank transform)         <base>__pilot  score-key zscore
  B  permutation z, weights permuted within hub-rank bins       <base>__pilot  score-key perm
  D  rank-based, housekeeping genes (HRT Atlas mouse) removed   <base>__nohk   score-key scores
Runs check_linked.py and check_anatomy.py --families ontology for every arm, then writes
results/pilot_stage1.tsv and prints the comparison. Success measures fixed in advance:
  * tissues whose own linked node ranks 1st (and top 3) of 74
  * Spearman between a linked node's link weight and its rank in its own tissue (currently -0.85; goal ~0)
  * ontology-family verdicts
Usage: pilot_stage1.py [--run-checks]
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
ARMS = [("A rank", "__pilot", "scores"), ("Z z-score", "__pilot", "zscore"), ("B perm z", "__pilot", "perm"), ("D no-HK", "__nohk", "scores")]


def run(cmd):
    r = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
    return (r.stdout + r.stderr).strip().splitlines()[-1:] if r.returncode else ["ok"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--run-checks", action="store_true"); a = ap.parse_args()
    pilot = json.load(open(ROOT / "data/_catalog/pilot_stage1.json"))
    bases = [e for e in pilot if e["arm"] != "D"]
    if a.run_checks:
        jobs = []
        for e in bases:
            for arm, suf, key in ARMS:
                name = e["base"] + suf; npz = ROOT / "psev_out" / f"{name}_GC_ranks.npz"
                if not npz.exists():
                    continue
                common = ["--npz", str(npz), "--tissue", e["family"], "--label", e["tissue"], "--name", name, "--score-key", key]
                jobs.append([str(HERE / "check_linked.py")] + common)
                jobs.append([str(HERE / "check_anatomy.py")] + common + ["--corrected", "--families", "ontology"])
        with ThreadPoolExecutor(6) as ex:
            bad = [j for j, r in zip(jobs, ex.map(run, jobs)) if r != ["ok"]]
        print(f"{len(jobs)} checks run, {len(bad)} failed", *bad[:3], sep="\n")
    # link weight of each linked node
    A = np.load(ROOT / "psev_out/_anatomy_psev.npz"); V = A["values"]; link = (V > 1e-4).mean(axis=0) > 0.5
    ni = pd.read_csv(ROOT.parent / "SPOKE_NASA/spoke_v_2/node_info.tsv", sep="\t"); an = ni[ni.Node_Type == "Anatomy"].reset_index(drop=True)
    lw = pd.Series((-V.mean(axis=0)[link]).argsort().argsort() + 1, index=an.Node_Name.astype(str).values[link])
    rows = []
    for e in bases:
        for arm, suf, key in ARMS:
            name = e["base"] + suf
            ldir = ROOT / "datasets" / name / ("corrected_linked" if key == "scores" else f"corrected_{key}_linked")
            odir = ROOT / "datasets" / name / ("corrected_ontology" if key == "scores" else f"corrected_{key}_ontology")
            if not (ldir / "summary.json").exists():
                continue
            L = json.load(open(ldir / "summary.json")); O = json.load(open(odir / "summary.json")) if (odir / "summary.json").exists() else {}
            tab = pd.read_csv(ldir / "linked_nodes.tsv", sep="\t") if (ldir / "linked_nodes.tsv").exists() else None
            own = None
            if tab is not None and not L.get("proxy") and tab.target.any():
                t = tab[tab.target].sort_values("pooled_rank").iloc[0]; own = (t.Node_Name, int(t.pooled_rank), int(lw.get(t.Node_Name, -1)))
            rows.append(dict(tissue=e["tissue"], arm=arm, linked_verdict=L["verdict"] + (" (proxy)" if L.get("proxy") else ""),
                             best_target_rank=L.get("best_target_rank"), own_node=own[0] if own else None,
                             own_node_rank=own[1] if own else None, own_node_link_weight=own[2] if own else None,
                             top_linked=", ".join(L.get("top_linked_nodes", [])[:3]),
                             onto_verdict=O.get("verdict"), onto_family_median_pct=round(O["family_median_pct"], 1) if O else None,
                             onto_family_rank=O.get("family_rank_among_organ_families")))
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("no pilot results yet")
    df.to_csv(ROOT / "results/pilot_stage1.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    for arm, _, _ in ARMS:
        d = df[df.arm == arm]
        if d.empty: continue
        own = d.dropna(subset=["own_node_rank"])
        rho = spearmanr(own.own_node_link_weight, own.own_node_rank).correlation if len(own) > 3 else float("nan")
        print(f"\n== {arm}: own linked node ranked 1st in {int((d.best_target_rank == 1).sum())}/{len(d)} tissues, top 3 in "
              f"{int((d.best_target_rank <= 3).sum())}/{len(d)}; Spearman(link weight, own rank) = {rho:+.2f}; "
              f"ontology verdicts {d.onto_verdict.value_counts().to_dict()}")
        print(d[["tissue", "linked_verdict", "best_target_rank", "own_node", "own_node_rank", "top_linked", "onto_verdict", "onto_family_median_pct", "onto_family_rank"]].to_string(index=False, max_colwidth=46))


if __name__ == "__main__":
    main()

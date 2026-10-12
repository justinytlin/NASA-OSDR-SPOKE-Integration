#!/usr/bin/env python3
"""Stage 2: the same ten tissues on other SPOKE builds (2025 full graph and its edge-type variants).

For each build tag (e.g. 2025, 2025noAeG, 2025upAuG) this runs, per pilot dataset:
  check_linked.py   --spoke <build> --anatomy-psev psev_out/_anatomy_psev_<tag>.npz --tag <tag>
  check_anatomy.py  --corrected --families ontology --spoke <build>   (output subfolder corrected_ontology)
on the npz written by gc_counts_to_psev.py --zip <build> --spoke <build> (names <base>__<tag>), then
tabulates own-node rank among that build's gene-linked Anatomy nodes next to the 2019 result
(datasets/<base>__pilot/corrected_linked) -> results/pilot_stage2.tsv.

Usage: pilot_stage2.py --tag 2025 --build /Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03 [--run-checks]
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent


def run(cmd):
    r = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
    return "ok" if r.returncode == 0 else (r.stdout + r.stderr).strip().splitlines()[-1]


def read(dir_):
    sj = dir_ / "summary.json"
    return json.load(open(sj)) if sj.exists() else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--build", required=True)
    ap.add_argument("--run-checks", action="store_true"); a = ap.parse_args()
    pilot = [e for e in json.load(open(ROOT / "data/_catalog/pilot_stage1.json")) if e["arm"] != "D"]
    apsev = ROOT / "psev_out" / f"_anatomy_psev_{a.tag}.npz"
    if a.run_checks:
        jobs = []
        for e in pilot:
            name = f"{e['base']}__{a.tag}"; npz = ROOT / "psev_out" / f"{name}_GC_ranks.npz"
            if not npz.exists():
                continue
            common = ["--npz", str(npz), "--tissue", e["family"], "--label", e["tissue"], "--name", name, "--spoke", a.build]
            jobs.append([str(HERE / "check_linked.py")] + common + ["--anatomy-psev", str(apsev), "--tag", a.tag])
            jobs.append([str(HERE / "check_anatomy.py")] + common + ["--corrected", "--families", "ontology"])
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(run, jobs))
        print(f"{len(jobs)} checks, {sum(r != 'ok' for r in res)} failed", *[r for r in res if r != "ok"][:4], sep="\n")
    # how many linked nodes does this build have, and which organ nodes are linked?
    A = np.load(apsev); V = A["values"]; link = (V > 1e-4).mean(axis=0) > 0.5
    ni = pd.read_csv(Path(a.build) / "node_info.tsv", sep="\t"); an = ni[ni.Node_Type == "Anatomy"].reset_index(drop=True)
    print(f"\n{a.tag}: {int(link.sum())} of {len(an)} Anatomy nodes are gene-linked (threshold: PSEV > 1e-4 for > 50 % of genes)")
    frac = (V > 1e-4).mean(axis=0); o = np.argsort(-frac)
    print("   fraction of genes above 1e-4 at the top nodes:", [(an.Node_Name[i][:18], round(float(frac[i]), 2)) for i in o[:6]], "...",
          [(an.Node_Name[i][:18], round(float(frac[i]), 2)) for i in o[int(link.sum()) - 2: int(link.sum()) + 3]])
    rows = []
    for e in pilot:
        base = e["base"]
        L19 = read(ROOT / "datasets" / f"{base}__pilot" / "corrected_linked"); O19 = read(ROOT / "datasets" / f"{base}__pilot" / "corrected_ontology")
        Lb = read(ROOT / "datasets" / f"{base}__{a.tag}" / f"corrected_linked_{a.tag}"); Ob = read(ROOT / "datasets" / f"{base}__{a.tag}" / "corrected_ontology")
        rows.append(dict(tissue=e["tissue"],
                         rank_2019=L19.get("best_target_rank") if L19 else None, n_linked_2019=L19.get("n_linked") if L19 else None,
                         verdict_2019=(L19["verdict"] + (" (proxy)" if L19.get("proxy") else "")) if L19 else None,
                         onto_2019=O19["verdict"] if O19 else None, onto_pct_2019=round(O19["family_median_pct"], 1) if O19 else None,
                         **{f"rank_{a.tag}": Lb.get("best_target_rank") if Lb else None, f"n_linked_{a.tag}": Lb.get("n_linked") if Lb else None,
                            f"verdict_{a.tag}": (Lb["verdict"] + (" (proxy)" if Lb.get("proxy") else "")) if Lb else None,
                            f"targets_{a.tag}": "; ".join(Lb.get("target_nodes", [])[:2]) if Lb else None,
                            f"top_{a.tag}": ", ".join(Lb.get("top_linked_nodes", [])[:3]) if Lb else None,
                            f"onto_{a.tag}": Ob["verdict"] if Ob else None, f"onto_pct_{a.tag}": round(Ob["family_median_pct"], 1) if Ob else None}))
    df = pd.DataFrame(rows)
    out = ROOT / "results" / f"pilot_stage2_{a.tag}.tsv"; df.to_csv(out, sep="\t", index=False)
    pd.set_option("display.width", 260); print(df.to_string(index=False, max_colwidth=40))


if __name__ == "__main__":
    main()

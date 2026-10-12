#!/usr/bin/env python3
"""Run the linked-node and ontology-family checks for every entry of a build manifest (4 in parallel).

Usage: check_build.py --tag 2025upAuG --build /Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only [--names a,b] [--force]
Needs psev_out/_anatomy_psev_<tag>.npz (extract_anatomy_psev.py) and psev_out/<name>_GC_ranks.npz.
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent


def run(cmd):
    r = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
    return "ok" if r.returncode == 0 else f"{cmd[1]}: " + (r.stdout + r.stderr).strip().splitlines()[-1]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--build", required=True)
    ap.add_argument("--names", default=None); ap.add_argument("--force", action="store_true"); a = ap.parse_args()
    want = set(a.names.split(",")) if a.names else None
    apsev = ROOT / "psev_out" / f"_anatomy_psev_{a.tag}.npz"; jobs = []
    for e in json.load(open(ROOT / "data/_catalog" / f"manifest_{a.tag}.json")):
        if want and e["name"] not in want: continue
        npz = ROOT / "psev_out" / f"{e['name']}_GC_ranks.npz"
        if not npz.exists(): continue
        common = ["--npz", str(npz), "--tissue", e["family"], "--label", e["tissue"], "--name", e["name"], "--spoke", a.build]
        if a.force or not (ROOT / "datasets" / e["name"] / f"corrected_linked_{a.tag}" / "summary.json").exists():
            jobs.append([str(HERE / "check_linked.py")] + common + ["--anatomy-psev", str(apsev), "--tag", a.tag])
        if a.force or not (ROOT / "datasets" / e["name"] / "corrected_ontology" / "summary.json").exists():
            jobs.append([str(HERE / "check_anatomy.py")] + common + ["--corrected", "--families", "ontology"])
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(run, jobs))
    bad = [r for r in res if r != "ok"]
    print(f"{len(jobs)} checks, {len(bad)} failed", *bad[:5], sep="\n")


if __name__ == "__main__":
    main()

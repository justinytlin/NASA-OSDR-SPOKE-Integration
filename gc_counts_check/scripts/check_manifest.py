#!/usr/bin/env python3
"""Run a check over many datasets (6 in parallel).

  (default)            raw check for manifest entries with a rank table and no summary yet
  --corrected          baseline-corrected, name families      -> datasets/<name>/baseline_corrected/
  --variant ontology   baseline-corrected, UBERON families    -> datasets/<name>/corrected_ontology/
  --variant linked     tissue rank among gene-linked nodes     -> datasets/<name>/corrected_linked/
The corrected variants cover every dataset folder (manifest entries and the first 14) whose
.npz carries scores and the uniform-weight profile.

Usage: check_manifest.py [--names a,b] [--force] [--corrected | --variant ontology|linked]
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def run(job):
    script, cmd, name = job
    r = subprocess.run([sys.executable, str(HERE / script)] + cmd, capture_output=True, text=True)
    line = next((l for l in r.stdout.splitlines() if l.startswith(name)), None)
    return line or f"{name}: FAILED {(r.stdout + r.stderr).strip().splitlines()[-1:]}"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--names", default=None); ap.add_argument("--force", action="store_true")
    ap.add_argument("--corrected", action="store_true"); ap.add_argument("--variant", choices=["ontology", "linked"], default=None)
    a = ap.parse_args()
    want = set(a.names.split(",")) if a.names else None
    jobs = []
    if a.corrected or a.variant:
        sub = {"ontology": "corrected_ontology", "linked": "corrected_linked"}.get(a.variant, "baseline_corrected")
        for sj in sorted((ROOT / "datasets").glob("*/summary.json")):
            s = json.load(open(sj)); name = s["name"]
            npz = ROOT / "psev_out" / f"{name}_GC_ranks.npz"
            if (want and name not in want) or not npz.exists() or ((sj.parent / sub / "summary.json").exists() and not a.force):
                continue
            base = ["--npz", str(npz), "--tissue", s["tissue"], "--label", s.get("tissue_label", s["tissue"]), "--name", name]
            if a.variant == "linked":
                jobs.append(("check_linked.py", base, name))
            else:
                jobs.append(("check_anatomy.py", base + ["--corrected"] + (["--families", "ontology"] if a.variant else []), name))
    else:
        for e in json.load(open(ROOT / "data/_catalog/manifest.json")):
            ranks = ROOT / "psev_out" / f"{e['name']}_GC_anatomy_ranks.tsv"
            if (want and e["name"] not in want) or not ranks.exists() or ((ROOT / "datasets" / e["name"] / "summary.json").exists() and not a.force):
                continue
            jobs.append(("check_anatomy.py", ["--ranks", str(ranks), "--tissue", e["family"], "--label", e["tissue"], "--name", e["name"]], e["name"]))
    with ThreadPoolExecutor(6) as ex:
        for line in ex.map(run, jobs):
            print(line, flush=True)


if __name__ == "__main__":
    main()

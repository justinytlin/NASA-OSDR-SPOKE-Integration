#!/usr/bin/env python3
"""Smoke test: run psev_pipeline + welch_meta against a small synthetic PSEV zip.

Builds a fake gene_psev.zip (20 groups, ~9 genes each, real node universe,
random PSEV values), runs the pipeline for all six studies, then welch_meta.
Outputs go to a scratch copy of the V2 dir so real runs are untouched.
"""

import io
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
V2 = BASE / "SPOKE_NASA/GeneLab_for_SPOKE/V2"
SCRATCH = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "synthetic_test"

sys.path.insert(0, str(HERE))
from psev_pipeline import (get_mapped_counts_and_diff_exp_dfs,
                           filter_same_direction,
                           load_homologene_mouse_to_human)

ACCESSIONS = ["GLDS-4", "GLDS-244", "GLDS-245", "GLDS-246", "GLDS-288", "GLDS-289"]

SCRATCH.mkdir(exist_ok=True)
work_v2 = SCRATCH / "V2"
if not work_v2.exists():
    work_v2.mkdir()
    for f in V2.glob("*differential_expression*.tsv"):
        (work_v2 / f.name).symlink_to(f)

node_info = pd.read_csv(BASE / "SPOKE_NASA/spoke_v_2/node_info.tsv", sep="\t")
spoke_genes = set(node_info[node_info.Node_Type == "Gene"].Node.values)
n_nodes = len(node_info)
gene_node_pos = {n: i for i, n in enumerate(node_info.Node.values)}

m2h = load_homologene_mouse_to_human(BASE / "psev_data/homologene_build68.data")

picked = []
for acc in ACCESSIONS:
    diff = next(f for f in V2.iterdir()
                if "differential_expression" in f.name and acc + "_" in f.name)
    exp = get_mapped_counts_and_diff_exp_dfs(diff, m2h, spoke_genes)
    samples = np.array([c for c in exp.columns if "Log2fc_" in c and c != "max_fc"])
    exp = filter_same_direction(exp, samples)
    kept = exp[exp.Space_over_Ground_Basal_same_sign == True].Node.values
    picked += list(np.random.default_rng(0).choice(kept, 30, replace=False))
genes = sorted(set(picked))
print(f"synthetic universe: {len(genes)} genes")

rng = np.random.default_rng(42)
fake_zip = SCRATCH / "gene_psev.zip"
per_group = int(np.ceil(len(genes) / 20))
with zipfile.ZipFile(fake_zip, "w", zipfile.ZIP_STORED) as zf:
    for g in range(20):
        chunk = genes[g * per_group:(g + 1) * per_group]
        tsv = "node_2_index\n" + "\n".join(str(gene_node_pos[n]) for n in chunk) + "\n"
        zf.writestr(f"gene_psevs/gene_group_{g}.tsv", tsv)
        mat = rng.random((len(chunk), n_nodes))
        buf = io.BytesIO()
        np.save(buf, mat)
        zf.writestr(f"gene_psevs/raw_psev_0_1_gene_group_{g}_sparse.npy", buf.getvalue())

for acc in ACCESSIONS:
    r = subprocess.run(
        [sys.executable, str(HERE / "psev_pipeline.py"), "--accession", acc,
         "--input", str(work_v2), "--zip", str(fake_zip)],
        capture_output=True, text=True)
    print(r.stdout.strip().splitlines()[-1] if r.returncode == 0 else r.stderr[-2000:])
    if r.returncode:
        sys.exit(1)

import os
env = dict(os.environ, WELCH_V2=str(work_v2))
r = subprocess.run(
    [sys.executable, str(HERE / "welch_meta.py"), "--out", str(SCRATCH / "meta_out")],
    capture_output=True, text=True, env=env)
print(r.stdout[-3000:] if r.returncode == 0 else r.stderr[-3000:])

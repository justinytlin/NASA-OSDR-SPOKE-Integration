#!/usr/bin/env python3
"""Production run: PSEV ranks for all six GLDS studies in one sweep.

Extracts each of the 20 PSEV gene-group matrices from gene_psev.zip once per
pass (2 passes total) and processes every study against it, instead of
re-extracting per study. Peak disk: zip + one 3.1 GB group file.

Usage: run_all.py [--zip PATH] [--input DIR] [--limit-groups N (debug)]
"""

import argparse
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent

sys.path.insert(0, str(HERE))
from psev_pipeline import (get_mapped_counts_and_diff_exp_dfs,
                           filter_same_direction,
                           load_homologene_mouse_to_human,
                           load_gene_indices_df, member_name,
                           get_order_rank_vector, get_rank_by_type_vector)

ACCESSIONS = ["GLDS-4", "GLDS-244", "GLDS-245", "GLDS-246", "GLDS-288", "GLDS-289"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(BASE / "SPOKE_NASA/GeneLab_for_SPOKE/V2"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    ap.add_argument("--limit-groups", type=int, default=None)
    args = ap.parse_args()

    input_dir = Path(args.input)
    save_str = "_".join(str(args.b).split("."))

    node_info_df = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    node_list = node_info_df.Node.values
    n_nodes = len(node_list)
    spoke_genes = set(node_info_df[node_info_df.Node_Type == "Gene"].Node.values)
    mouse_to_human = load_homologene_mouse_to_human(args.homologene)

    # ---- prep every study ------------------------------------------------
    studies = {}
    for acc in ACCESSIONS:
        diff_file = [f for f in input_dir.iterdir()
                     if "differential_expression" in f.name and acc + "_" in f.name][0]
        exp_df = get_mapped_counts_and_diff_exp_dfs(diff_file, mouse_to_human, spoke_genes)
        samples = np.array([c for c in exp_df.columns if "Log2fc_" in c and c != "max_fc"])
        exp_df = filter_same_direction(exp_df, samples)
        kept = exp_df[exp_df.Space_over_Ground_Basal_same_sign == True]
        studies[acc] = {"kept": kept, "samples": samples}
        print(f"[{acc}] genes kept: {len(kept)}, comparisons: {len(samples)}", flush=True)

    with zipfile.ZipFile(args.zip) as zf:
        gene_indices_df = load_gene_indices_df(zf, node_list)
        rounds = sorted(gene_indices_df.Round.unique())
        if args.limit_groups:
            rounds = rounds[:args.limit_groups]

        for acc, st in studies.items():
            gi = gene_indices_df.merge(st["kept"][["Node"]], on="Node")
            st["gi"] = gi
            st["total"] = len(gi)
            st["exp"] = st["kept"].merge(gi[["Node", "Round", "round_index"]], on="Node")
            st["sum"] = np.zeros(n_nodes)
            st["sumsq"] = np.zeros(n_nodes)
            print(f"[{acc}] genes with PSEVs: {st['total']}", flush=True)

        tmpdir = Path(tempfile.mkdtemp(prefix="psev_"))

        def load_group(r):
            name = member_name(zf, f"raw_psev_{save_str}_gene_group_{r}_sparse.npy")
            target = tmpdir / f"g{r}.npy"
            with zf.open(name) as src, open(target, "wb") as dst:
                while True:
                    chunk = src.read(1 << 24)
                    if not chunk:
                        break
                    dst.write(chunk)
            return target

        # ---- pass 1: per-study mean/std over its seen genes -------------
        for r in rounds:
            target = load_group(r)
            mat = np.load(target, mmap_mode="r")
            for acc, st in studies.items():
                idx = st["gi"][st["gi"].Round == r].round_index.values
                if not len(idx):
                    continue
                sub = np.asarray(mat[np.sort(idx)], dtype=np.float64)
                st["sum"] += sub.sum(axis=0)
                st["sumsq"] += (sub ** 2).sum(axis=0)
                del sub
            del mat
            target.unlink()
            print(f"pass1 group {r} done", flush=True)

        for acc, st in studies.items():
            avg = st["sum"] / float(st["total"])
            std = np.sqrt(st["sumsq"] / float(st["total"]) - avg ** 2)
            std[std == 0] = np.nan
            st["avg"], st["std"] = avg, std
            st["acc_psev"] = np.zeros((len(st["samples"]), n_nodes))

        # ---- pass 2: z-score, rank rows, dot with FC --------------------
        order = np.arange(n_nodes)
        for r in rounds:
            target = load_group(r)
            mat = np.load(target, mmap_mode="r")
            for acc, st in studies.items():
                gi_r = st["gi"][st["gi"].Round == r]
                idx = gi_r.round_index.values
                if not len(idx):
                    continue
                idx_sorted = np.sort(idx)
                sub = np.asarray(mat[idx_sorted], dtype=np.float64)
                z = np.nan_to_num((sub - st["avg"]) / st["std"])
                del sub
                ranked = np.empty_like(z)
                for i, row in enumerate(z):
                    ranked[i] = get_order_rank_vector(row, order)
                del z
                fc = (st["exp"][st["exp"].Round == r].set_index("round_index")
                      .loc[idx_sorted, st["samples"]].values)
                st["acc_psev"] += fc.T @ ranked
                del ranked
            del mat
            target.unlink()
            print(f"pass2 group {r} done", flush=True)

    # ---- final ranking + save -------------------------------------------
    order = np.arange(n_nodes)
    for acc, st in studies.items():
        sample_psev = st["acc_psev"] / float(st["total"])
        ranked_out = np.array([get_order_rank_vector(row, order) for row in sample_psev])
        result_df = pd.concat(
            (node_info_df, pd.DataFrame(ranked_out.T, columns=st["samples"])), axis=1)
        for sample in st["samples"]:
            result_df.loc[:, "Rank_by_type_%s" % sample] = get_rank_by_type_vector(
                result_df[sample].values, node_info_df.Node_Type.values)
        out = input_dir / f"{acc}_ranks_and_rank_by_type_for_meta{save_str}.tsv"
        result_df.to_csv(out, sep="\t", index=False)
        print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()

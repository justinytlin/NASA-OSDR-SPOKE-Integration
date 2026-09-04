#!/usr/bin/env python3
"""Replication of Nelson et al., Life 2021, 11:42 — GeneLab FC -> PSEV node ranks.

Faithful port of baranzini-lab/SPOKE_NASA genelab_fc_to_psev.py with three changes:
  * git merge conflict resolved (single-study -i/--accession form kept)
  * mouse->human orthologs from frozen NCBI HomoloGene build 68 (the 2020-era
    MGI report has since dropped the "HomoloGene ID" column)
  * PSEV gene-group matrices are streamed out of gene_psev.zip one at a time
    (peak disk ~zip + one group, instead of a 60 GB full extraction)

Usage:
  psev_pipeline.py --accession GLDS-4 [--prep-only]
      [--zip PSEV_ZIP] [--spoke SPOKE_DIR] [--input INPUT_DIR]
      [--homologene FILE] [-b 0.1]

Outputs <input>/<accession>_ranks_and_rank_by_type_for_meta<b>.tsv, identical in
shape to the original: node_info columns + per-comparison overall rank and
Rank_by_type_* columns.
"""

import argparse
import itertools
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def get_order_rank_vector(test_vector, ground_truth_order):
    temp = test_vector.argsort()
    ranks = np.empty(len(test_vector), int)
    ranks[temp] = np.arange(len(test_vector))
    return ranks[ground_truth_order]


def get_rank_by_type_vector(vector, type_list):
    rank_list = np.zeros(len(type_list))
    for node_type in set(type_list):
        node_type_index = np.arange(len(type_list))[type_list == node_type]
        rank_list[node_type_index] = len(node_type_index) - get_order_rank_vector(
            vector[node_type_index], np.arange(len(node_type_index)))  # low is top
    return rank_list


def load_homologene_mouse_to_human(path):
    """HomoloGene build68: HID, taxid, entrez, symbol, ... -> mouse/human Entrez pairs."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["HID", "taxid", "entrez", "symbol", "gi", "acc"],
                     usecols=["HID", "taxid", "entrez", "symbol"])
    mouse = df[df.taxid == 10090][["HID", "entrez", "symbol"]].rename(
        columns={"entrez": "ENTREZID", "symbol": "Mouse_Symbol"})
    human = df[df.taxid == 9606][["HID", "entrez"]].rename(
        columns={"entrez": "Human_EntrezGene ID"})
    return pd.merge(mouse, human, on="HID")


class PsevStore:
    """Gene-PSEV store: either the Zenodo gene_psev.zip or a directory laid out
    the same way (psev_build/make_gene_psevs.py output). Members are located by
    filename suffix, so 'gene_psevs/' prefixes and __MACOSX junk are ignored."""

    def __init__(self, path):
        self.path = Path(path)
        self.zf = None
        if self.path.is_dir():
            self.names = [str(f.relative_to(self.path))
                          for f in self.path.rglob("*") if f.is_file()]
        else:
            self.zf = zipfile.ZipFile(self.path)
            self.names = self.zf.namelist()
        # ignore macOS AppleDouble sidecars (._foo) that exFAT/zip archives carry
        self.names = [n for n in self.names if not Path(n).name.startswith("._")]

    def close(self):
        if self.zf:
            self.zf.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def member(self, suffix):
        for n in self.names:
            if n.endswith(suffix) and "__MACOSX" not in n:
                return n
        raise KeyError(f"{suffix} not found in {self.path}")

    def open(self, suffix):
        name = self.member(suffix)
        return self.zf.open(name) if self.zf else open(self.path / name, "rb")

    def groups(self):
        gs = sorted(int(re.search(r"gene_group_(\d+)\.tsv$", n).group(1))
                    for n in self.names if re.search(r"gene_group_(\d+)\.tsv$", n)
                    and "__MACOSX" not in n)
        return gs

    def load_group(self, r, save_str):
        """Dense (genes_in_group, n_nodes) matrix for one group."""
        name = self.member(f"raw_psev_{save_str}_gene_group_{r}_sparse.npy")
        if not self.zf:
            return np.load(self.path / name, allow_pickle=False, mmap_mode="r")
        tmpdir = Path(tempfile.mkdtemp(prefix="psev_"))
        target = tmpdir / f"g{r}.npy"
        with self.zf.open(name) as src, open(target, "wb") as dst:
            while True:
                chunk = src.read(1 << 24)
                if not chunk:
                    break
                dst.write(chunk)
        mat = np.load(target, allow_pickle=False, mmap_mode=None)
        target.unlink()
        return mat


def load_gene_indices_df(store, node_list):
    frames = []
    for group_index in store.groups():
        with store.open(f"gene_group_{group_index}.tsv") as fh:
            df = pd.read_csv(fh, sep="\t", header=0, index_col=False)
        df.loc[:, "Round"] = group_index
        df.loc[:, "round_index"] = np.arange(len(df))
        frames.append(df)
    gene_indices_df = pd.concat(frames, axis=0)
    gene_indices_df.loc[:, "Node"] = [node_list[i] for i in gene_indices_df.node_2_index.values]
    return gene_indices_df


def get_mapped_counts_and_diff_exp_dfs(diff_path, mouse_to_human, spoke_genes):
    exp_df = pd.read_csv(diff_path, sep="\t", header=0, index_col=False, low_memory=False)
    exp_df = exp_df[exp_df.ENTREZID.notna()]
    # newer GLbulkRNAseq tables carry multi-mapped ids like "54204|100043580"
    exp_df.loc[:, "ENTREZID"] = exp_df.ENTREZID.astype(str).str.split("|")
    exp_df = exp_df.explode("ENTREZID")
    exp_df.loc[:, "ENTREZID"] = exp_df.ENTREZID.astype(float).astype(int)
    # map from ensembl+mouse-entrez to human entrez
    map_from_ensembl = pd.merge(
        exp_df[["ENSEMBL", "ENTREZID"]],
        mouse_to_human[["ENTREZID", "Human_EntrezGene ID"]], on="ENTREZID")
    map_from_ensembl.loc[:, "Human_EntrezGene ID"] = (
        map_from_ensembl["Human_EntrezGene ID"].values.astype(int).astype(str))
    map_from_ensembl = map_from_ensembl[
        map_from_ensembl["Human_EntrezGene ID"].isin(spoke_genes)
    ].rename(columns={"Human_EntrezGene ID": "Node"})
    exp_df = pd.merge(map_from_ensembl, exp_df, on=["ENSEMBL", "ENTREZID"])
    keep = [c for c in exp_df.columns
            if re.match(r"Node|P.value_|Adj.p.value_|Log2fc_|LRT.p.value", c)]
    exp_df = exp_df[keep].drop_duplicates()
    # mean for genes seen more than once; drop all-zero fc
    exp_df = exp_df.groupby("Node").mean(numeric_only=True).reset_index()
    fc_cols = [c for c in exp_df.columns if "Log2fc_" in c]
    exp_df.loc[:, "max_fc"] = np.max(exp_df[fc_cols].values, axis=1)
    exp_df = exp_df[exp_df.max_fc != 0]
    return exp_df


def check_sign(exp_df, groups, space_ground_basal, group_type_1, group_type_2):
    group_1_list = [groups[i] for i, t in enumerate(space_ground_basal) if group_type_1 == t]
    group_2_list = [groups[i] for i, t in enumerate(space_ground_basal) if group_type_2 == t]
    group_1_over_2 = ["Log2fc_(%s)v(%s)" % t for t in itertools.product(group_1_list, group_2_list)]
    group_1_over_2 = [c for c in group_1_over_2 if c in exp_df.columns]
    if not group_1_over_2:
        raise RuntimeError(f"no {group_type_1} over {group_type_2} columns found")
    vals = exp_df[group_1_over_2].values
    exp_df.loc[:, "%s_over_%s_same_sign" % (group_type_1, group_type_2)] = (
        (np.sum(vals > 0, axis=1) == len(group_1_over_2)) |
        (np.sum(vals < 0, axis=1) == len(group_1_over_2)))
    exp_df.loc[:, "%s_over_%s_pos" % (group_type_1, group_type_2)] = (
        np.sum(vals > 0, axis=1) == len(group_1_over_2))
    return exp_df


def filter_same_direction(exp_df, samples):
    groups = np.array(list(set(np.ravel([s[8:-1].split(")v(") for s in samples]))))
    space_ground_basal = np.array(
        [re.findall(r"SPACE|GROUND|BASAL", g, re.IGNORECASE)[0].capitalize() for g in groups])
    control_groups = np.setdiff1d(space_ground_basal, ["Space"])
    for control_group in control_groups:
        exp_df = check_sign(exp_df, groups, space_ground_basal, "Space", control_group)
    same = np.all(
        exp_df[["Space_over_%s_same_sign" % c for c in control_groups]].values == True, axis=1)
    pos = exp_df[["Space_over_%s_pos" % c for c in control_groups]].values
    exp_df.loc[:, "Space_over_Ground_Basal_same_sign"] = same & (
        np.all(pos == True, axis=1) | np.all(pos == False, axis=1))
    return exp_df


def compute_rank_table(kept, samples, node_info_df, zip_path, b=0.1, limit_groups=None):
    """Weighted-PSEV node ranks for one study (passes 1+2 of the pipeline).

    zip_path may be the Zenodo gene_psev.zip or a directory with the same
    layout (see psev_build/make_gene_psevs.py).
    kept: DataFrame with a 'Node' column (SPOKE gene ids, str) plus one numeric
    weight column per entry in `samples` (log2 fold changes). Returns
    node_info_df + per-sample overall-rank and Rank_by_type_* columns.
    """
    save_str = "_".join(str(b).split("."))
    node_list = node_info_df.Node.values
    n_nodes = len(node_list)
    samples = list(samples)

    with PsevStore(zip_path) as zf:
        gene_indices_df = load_gene_indices_df(zf, node_list)
        gene_indices_df = gene_indices_df.merge(kept[["Node"]], on="Node")
        total_genes = len(gene_indices_df)
        print(f"  genes with PSEVs: {total_genes}")

        exp_mat_df = kept.merge(
            gene_indices_df[["Node", "Round", "round_index"]], on="Node")

        def load_group(r):
            return zf.load_group(r, save_str)

        rounds = sorted(gene_indices_df.Round.unique())
        if limit_groups:
            rounds = rounds[:limit_groups]
        # pass 1: mean/std over seen genes
        s = np.zeros(n_nodes, dtype=np.float64)
        ss = np.zeros(n_nodes, dtype=np.float64)
        for r in rounds:
            idx = gene_indices_df[gene_indices_df.Round == r].round_index.values
            mat = load_group(r)[idx].astype(np.float64)
            s += mat.sum(axis=0)
            ss += (mat ** 2).sum(axis=0)
            del mat
            print(f"  pass1 group {r} done")
        avg_rank = s / float(total_genes)
        std_rank = np.sqrt(ss / float(total_genes) - avg_rank ** 2)
        std_rank[std_rank == 0] = np.nan

        # pass 2: z-score, rank rows, dot with FC
        sample_psev = np.zeros((len(samples), n_nodes), dtype=np.float64)
        for r in rounds:
            gi = gene_indices_df[gene_indices_df.Round == r]
            idx = gi.round_index.values
            z = np.nan_to_num((load_group(r)[idx].astype(np.float64) - avg_rank) / std_rank)
            ranked = np.array([get_order_rank_vector(row, np.arange(n_nodes)) for row in z],
                              dtype=np.float64)
            del z
            fc = (exp_mat_df[exp_mat_df.Round == r]
                  .set_index("round_index").loc[idx, samples].values)
            sample_psev += fc.T @ ranked
            del ranked
            print(f"  pass2 group {r} done")

        sample_psev /= float(total_genes)

    ranked_out = np.array([get_order_rank_vector(row, np.arange(n_nodes)) for row in sample_psev])
    result_df = pd.concat(
        (node_info_df, pd.DataFrame(ranked_out.T, columns=samples)), axis=1)
    for sample in samples:
        result_df.loc[:, "Rank_by_type_%s" % sample] = get_rank_by_type_vector(
            result_df[sample].values, node_info_df.Node_Type.values)
    return result_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True, help="e.g. GLDS-4")
    ap.add_argument("--input", default=str(BASE / "SPOKE_NASA/GeneLab_for_SPOKE/V2"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    ap.add_argument("--prep-only", action="store_true",
                    help="stop after gene mapping/filtering; print summary")
    args = ap.parse_args()

    input_dir = Path(args.input)
    save_str = "_".join(str(args.b).split("."))

    node_info_df = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    node_list = node_info_df.Node.values
    spoke_genes = set(node_info_df[node_info_df.Node_Type == "Gene"].Node.values)

    diff_file = [f for f in os.listdir(input_dir)
                 if "differential_expression" in f and args.accession + "_" in f
                 and not f.startswith("._")][0]
    print(f"[{args.accession}] input: {diff_file}")

    mouse_to_human = load_homologene_mouse_to_human(args.homologene)
    exp_df = get_mapped_counts_and_diff_exp_dfs(input_dir / diff_file, mouse_to_human, spoke_genes)
    samples = np.array([c for c in exp_df.columns if "Log2fc_" in c and c != "max_fc"])
    print(f"  mapped SPOKE genes: {len(exp_df)}; comparisons: {len(samples)}")

    exp_df = filter_same_direction(exp_df, samples)
    kept = exp_df[exp_df.Space_over_Ground_Basal_same_sign == True]
    print(f"  same-direction genes kept: {len(kept)}")
    if args.prep_only:
        return

    result_df = compute_rank_table(kept, samples, node_info_df, args.zip, b=args.b)

    out = input_dir / f"{args.accession}_ranks_and_rank_by_type_for_meta{save_str}.tsv"
    result_df.to_csv(out, sep="\t", index=False)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

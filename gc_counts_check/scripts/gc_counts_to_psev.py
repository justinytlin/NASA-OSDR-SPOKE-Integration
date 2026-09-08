#!/usr/bin/env python3
"""Ground-control expression -> per-sample SPOKE (2019 PSEV) node ranks, batched.

Same computation as replication/counts_to_psev.py (which wraps
replication/psev_pipeline.compute_rank_table), but all datasets are processed
in one streaming pass over the 25.6 GB Zenodo gene_psev.zip instead of one
pass per dataset. Per dataset the maths is identical to the Nelson 2021 port:

  pass 1: mean/std of each node's PSEV value over the dataset's genes
  pass 2: z-score each gene PSEV, rank the nodes within the gene, then
          sample_psev = weights.T @ ranks / n_genes
  final:  overall rank per sample + Rank_by_type_ per sample

Per-gene, per-sample weights:
  RNA-seq counts : log2(CPM + 1), genes with CPM >= 1 in >= 1 kept sample
  Affymetrix     : GeneLab RMA-normalized log2 probeset expression, collapsed
                   to one probeset per gene (highest mean; "MaxMean")
Gene ids: ENSEMBL -> mouse Entrez (GeneLab annotation) -> human Entrez
(HomoloGene build 68) -> SPOKE v2 Gene node.

OSD-4 (Affymetrix thymus) was excluded at the user's request; the 'array'
kind is kept for reference.

Usage:
  gc_counts_to_psev.py [--only OSD-564 ...] [--prep-only] [--threads 4]
  gc_counts_to_psev.py --counts PATH.csv --name OSD-xxx --tissue liver [--samples _GC_]   (any dataset)
Outputs (GC Counts Check/psev_out/):
  <name>_GC_weights.tsv                        prepared gene x sample weights
  <name>_GC_counts_ranks_and_rank_by_type_for_meta0_1.tsv   rank table
"""
import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                 # GC Counts Check/
BASE = ROOT.parent                 # SPOKE/
sys.path.insert(0, str(BASE / "replication"))
from psev_pipeline import (PsevStore, load_gene_indices_df, load_homologene_mouse_to_human,
                           get_order_rank_vector, get_rank_by_type_vector)
from counts_to_psev import load_ensembl_to_entrez

DATA = ROOT / "data"
OUT = ROOT / "psev_out"

DATASETS = [
    dict(name="OSD-564", tissue="hippocampus", kind="counts",
         file="GLDS-569_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-244", tissue="thymus", kind="counts",
         file="GLDS-244_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-245", tissue="liver", kind="counts",
         file="GLDS-245_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-246", tissue="spleen", kind="counts",
         file="GLDS-246_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-288", tissue="spleen", kind="counts",
         file="GLDS-288_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-289", tissue="thymus", kind="counts",
         file="GLDS-289_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    # RR-3 (SpaceX-8) payload
    dict(name="OSD-137", tissue="liver", kind="counts",
         file="GLDS-137_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-161", tissue="adrenal gland", kind="counts",
         file="GLDS-161_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-162", tissue="eye", kind="counts",
         file="GLDS-162_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-163", tissue="kidney", kind="counts",
         file="GLDS-163_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-168", tissue="liver", kind="counts",   # RR-1 + RR-3 study; RR-3 GC samples only
         file="GLDS-168_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_RR3_GC_"]),
    dict(name="OSD-194", tissue="retina", kind="counts",
         file="GLDS-194_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-270", tissue="heart", kind="counts",
         file="GLDS-270_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
    dict(name="OSD-352", tissue="brain", kind="counts",
         file="GLDS-352_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", samples=["_GC_"]),
]
RR3 = ["OSD-137", "OSD-161", "OSD-162", "OSD-163", "OSD-168", "OSD-194", "OSD-270", "OSD-352"]


def mouse_ensembl_to_spoke(e2e, m2h, spoke_genes):
    """ENSEMBL -> SPOKE human gene node table (columns ENSEMBL, Node)."""
    e2e = e2e.copy()
    e2e["ENTREZID"] = e2e.ENTREZID.astype(str).str.split("|")
    e2e = e2e.explode("ENTREZID")
    e2e = e2e[e2e.ENTREZID.str.fullmatch(r"\d+(\.0)?")]
    e2e["ENTREZID"] = e2e.ENTREZID.astype(float).astype(int)
    e2e = e2e.merge(m2h[["ENTREZID", "Human_EntrezGene ID"]], on="ENTREZID")
    e2e["Node"] = e2e["Human_EntrezGene ID"].astype(int).astype(str)
    return e2e[e2e.Node.isin(spoke_genes)][["ENSEMBL", "Node"]].drop_duplicates()


def prep_counts(path, filters, min_cpm=1.0):
    df = pd.read_csv(path, header=0, index_col=0)
    keep = [c for c in df.columns if any(s in c for s in filters)]
    df = df[keep]
    cpm = df / df.sum(axis=0) * 1e6
    expressed = (cpm >= min_cpm).any(axis=1)
    w = np.log2(cpm[expressed] + 1.0)
    w.index.name = "ENSEMBL"
    info = f"{df.shape[0]} genes, {len(keep)} GC samples, {int(expressed.sum())} with CPM>={min_cpm}"
    return w.reset_index(), info


def prep_array(path, filters):
    df = pd.read_csv(path, header=0, low_memory=False)
    samples = [c for c in df.columns if any(s in c for s in filters)]
    df = df[df.ENSEMBL.notna()]
    # MaxMean collapse: keep the probeset with the highest mean over kept samples
    df = df.assign(_m=df[samples].mean(axis=1)).sort_values("_m", ascending=False)
    df = df.drop_duplicates("ENSEMBL")
    w = df[["ENSEMBL"] + samples].reset_index(drop=True)
    info = f"{len(w)} genes (probesets collapsed MaxMean), {len(samples)} GC samples"
    return w, info


def prepare(ds, e2e_map, spoke_genes, m2h):
    path = Path(ds["file"]) if Path(str(ds["file"])).is_absolute() else DATA / ds["file"]
    if ds["kind"] == "counts":
        w, info = prep_counts(path, ds["samples"])
        emap = e2e_map
    else:
        w, info = prep_array(path, ds["samples"])
        ann = pd.read_csv(path, usecols=["ENSEMBL", "ENTREZID"], low_memory=False)
        emap = mouse_ensembl_to_spoke(ann[ann.ENTREZID.notna()], m2h, spoke_genes)
    samples = [c for c in w.columns if c != "ENSEMBL"]
    kept = w.merge(emap, on="ENSEMBL").groupby("Node")[samples].mean().reset_index()
    print(f"[{ds['name']} {ds['tissue']}] {ds['file']}: {info}; mapped SPOKE genes: {len(kept)}",
          flush=True)
    return kept, samples


def rank_rows(z):
    """rank of each element within its row (0 = smallest), matches get_order_rank_vector."""
    order = z.argsort(axis=1)
    ranks = np.empty(z.shape, dtype=np.float64)
    np.put_along_axis(ranks, order, np.broadcast_to(np.arange(z.shape[1], dtype=np.float64), z.shape), axis=1)
    return ranks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=None, help="dataset name(s); 'RR3' expands to the RR-3 set")
    ap.add_argument("--counts", default=None, help="ad hoc dataset: path to a GeneLab unnormalized counts CSV")
    ap.add_argument("--name", default=None, help="ad hoc dataset: label used in output file names")
    ap.add_argument("--tissue", default="unknown", help="ad hoc dataset: tissue label (informational)")
    ap.add_argument("--samples", action="append", default=None,
                    help="ad hoc dataset: substring(s) selecting sample columns (default _GC_)")
    ap.add_argument("--prep-only", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=128)
    ap.add_argument("--limit-groups", type=int, default=None)
    ap.add_argument("--zip", default=str(BASE / "psev_data/gene_psev.zip"))
    ap.add_argument("--spoke", default=str(BASE / "SPOKE_NASA/spoke_v_2"))
    ap.add_argument("--gene-map", default=str(BASE / "psev_data/OSD-564_rna_seq_differential_expression.csv"))
    ap.add_argument("--homologene", default=str(BASE / "psev_data/homologene_build68.data"))
    ap.add_argument("-b", type=float, default=0.1)
    args = ap.parse_args()
    save_str = "_".join(str(args.b).split("."))
    OUT.mkdir(exist_ok=True)

    node_info = pd.read_csv(Path(args.spoke) / "node_info.tsv", sep="\t")
    node_list = node_info.Node.values
    n_nodes = len(node_list)
    spoke_genes = set(node_info[node_info.Node_Type == "Gene"].Node.values)
    m2h = load_homologene_mouse_to_human(args.homologene)
    e2e_map = mouse_ensembl_to_spoke(load_ensembl_to_entrez(args.gene_map), m2h, spoke_genes)

    if args.counts:
        cp = Path(args.counts)
        todo = [dict(name=args.name or cp.stem.split("_")[0], tissue=args.tissue, kind="counts",
                     file=cp.resolve(), samples=args.samples or ["_GC_"])]
    else:
        only = None
        if args.only:
            only = [n for o in args.only for n in (RR3 if o.upper() == "RR3" else [o])]
        todo = [d for d in DATASETS if not only or d["name"] in only]
    prepared = {}
    for ds in todo:
        kept, samples = prepare(ds, e2e_map, spoke_genes, m2h)
        kept.to_csv(OUT / f"{ds['name']}_GC_weights.tsv", sep="\t", index=False)
        prepared[ds["name"]] = (kept, samples)
    if args.prep_only:
        return

    t0 = time.time()
    with PsevStore(args.zip) as zf, ThreadPoolExecutor(max_workers=args.threads) as pool:
        gi_all = load_gene_indices_df(zf, node_list)
        rounds = sorted(gi_all.Round.unique())
        if args.limit_groups:
            rounds = rounds[:args.limit_groups]
        st = {}
        for name, (kept, samples) in prepared.items():
            gi = gi_all.merge(kept[["Node"]], on="Node")
            exp = kept.merge(gi[["Node", "Round", "round_index"]], on="Node")
            st[name] = dict(gi=gi, exp=exp, samples=samples, total=len(gi),
                            s=np.zeros(n_nodes), ss=np.zeros(n_nodes),
                            psev=np.zeros((len(samples), n_nodes)))
            print(f"  {name}: genes with PSEVs: {len(gi)}", flush=True)

        # ---- pass 1: per-dataset mean / std over its genes ----
        for r in rounds:
            mat = zf.load_group(r, save_str)
            for name, d in st.items():
                idx = d["gi"][d["gi"].Round == r].round_index.values
                for i in range(0, len(idx), args.chunk):
                    sub = mat[idx[i:i + args.chunk]].astype(np.float64)
                    d["s"] += sub.sum(axis=0)
                    d["ss"] += (sub ** 2).sum(axis=0)
            del mat
            print(f"  pass1 group {r} done ({time.time() - t0:.0f}s)", flush=True)
        for d in st.values():
            d["avg"] = d["s"] / float(d["total"])
            with np.errstate(invalid="ignore"):
                d["std"] = np.sqrt(d["ss"] / float(d["total"]) - d["avg"] ** 2)
            d["std"][d["std"] == 0] = np.nan

        # ---- pass 2: z-score, rank within gene, weight by expression ----
        def work(mat, idx, d):
            z = np.nan_to_num((mat[idx].astype(np.float64) - d["avg"]) / d["std"])
            ranked = rank_rows(z)
            del z
            fc = d["exp"][d["exp"].Round == r].set_index("round_index").loc[idx, d["samples"]].values
            return fc.T @ ranked

        for r in rounds:
            mat = zf.load_group(r, save_str)
            futs = []
            for name, d in st.items():
                idx = d["gi"][d["gi"].Round == r].round_index.values
                for i in range(0, len(idx), args.chunk):
                    futs.append((name, pool.submit(work, mat, idx[i:i + args.chunk], d)))
            for name, f in futs:
                st[name]["psev"] += f.result()
            del mat
            print(f"  pass2 group {r} done ({time.time() - t0:.0f}s)", flush=True)

    for name, d in st.items():
        sp = d["psev"] / float(d["total"])
        ranked_out = np.array([get_order_rank_vector(row, np.arange(n_nodes)) for row in sp])
        res = pd.concat((node_info, pd.DataFrame(ranked_out.T, columns=d["samples"])), axis=1)
        for s in d["samples"]:
            res["Rank_by_type_" + s] = get_rank_by_type_vector(res[s].values, node_info.Node_Type.values)
        out = OUT / f"{name}_GC_counts_ranks_and_rank_by_type_for_meta{save_str}.tsv"
        res.to_csv(out, sep="\t", index=False)
        print(f"  wrote {out}", flush=True)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

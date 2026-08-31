#!/usr/bin/env python3
"""Project a NASA OSDR spaceflight gene signature onto the SPOKE knowledge graph.

Dataset: OSD-137 (Rodent Research-3 CASIS) — mouse liver RNA-seq,
Space Flight vs Ground Control, from OSDR's processed visualization table:
  https://visualization.osdr.nasa.gov/biodata/api/v2/query/data/
      ?id.accession=OSD-137&file.datatype=visualization table&format=csv

Pipeline:
  1. select spaceflight-responsive candidate genes (nominal p < 0.01)
  2. map mouse -> human orthologs (MGI homology table via --orthologs,
     falling back to the uppercase-symbol heuristic without it)
  3. look each gene up in SPOKE; pull its Disease + Pathway neighborhood
  4. aggregate: which diseases/pathways the spaceflight signature lands on

Outputs (into ./osdr_projection/):
  gene_matches.csv   per-gene: mouse stats + whether SPOKE knows the gene
  diseases.csv       diseases ranked by number of signature genes hitting them
  pathways.csv       pathways ranked the same way
  edges.csv          every gene->disease / gene->pathway edge (for graph viz)

Usage:
  python3 examples/osdr_spoke_projection.py TABLE.csv [--out DIR] [--p 0.01]
      [--fdr 0.05] [--top N]

  --p     select by nominal p-value (default mode, cutoff 0.01)
  --fdr   select by adjusted p instead (e.g. --fdr 0.05)
  --top   cap the gene list at the N most significant (default: no cap)
  --out   output directory (default: osdr_projection/ next to spoke_api/)
  --orthologs FILE   MGI HOM_MouseHumanSequence.rpt for real ortholog mapping
  --human            input genes are already human symbols (skip mapping)
  --contrast NAME    DE contrast, default "(Space Flight)v(Ground Control)"
"""

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spoke_api import SpokeClient, split_graph
from spoke_api.orthologs import load_mgi, map_gene

DEFAULT_CONTRAST = "(Space Flight)v(Ground Control)"


def load_signature(table_path, p_cutoff, fdr_cutoff, top, contrast):
    """DE candidates with a real gene symbol, by nominal p or FDR."""
    r = csv.reader(open(table_path))
    hdr = [h.split("/")[-1] for h in next(r)]
    ix = {h: i for i, h in enumerate(hdr)}
    sym, lfc = ix["SYMBOL"], ix[f"Log2fc_{contrast}"]
    ent = ix.get("ENTREZID")
    p, adj = ix[f"P.value_{contrast}"], ix[f"Adj.p.value_{contrast}"]
    genes = []
    for row in r:
        try:
            pv, av, lv = float(row[p]), float(row[adj]), float(row[lfc])
        except ValueError:
            continue
        s = row[sym]
        if not s or s == "NaN":
            continue
        if (fdr_cutoff is not None and av < fdr_cutoff) or \
           (fdr_cutoff is None and pv < p_cutoff):
            entrez = row[ent] if ent is not None else ""
            if entrez and entrez != "NaN":
                entrez = entrez.split(".")[0]  # "12345.0" -> "12345"
            else:
                entrez = ""
            genes.append({"mouse": s, "mouse_entrez": entrez,
                          "log2fc": lv, "p": pv, "adj_p": av})
    genes.sort(key=lambda g: g["p"])
    return genes[:top] if top else genes


def add_orthologs(genes, mgi_path, human=False):
    """Attach human ortholog symbols; MGI table, identity for human data, else uppercase."""
    if human:
        for g in genes:
            g["human"] = g["mouse"]
            g["homology"] = "human"
        return genes
    if not mgi_path:
        for g in genes:
            g["human"] = g["mouse"].upper()
            g["homology"] = "uppercase-heuristic"
        return genes
    by_entrez, by_symbol = load_mgi(mgi_path)
    out = []
    for g in genes:
        hits = map_gene(g["mouse_entrez"], g["mouse"], by_entrez, by_symbol)
        if not hits:
            g["human"], g["homology"] = "", "no-ortholog"
            out.append(g)
            continue
        for h in hits:
            gg = dict(g)
            gg["human"] = h["symbol"]
            gg["homology"] = f"{h['n_mouse']}:{h['n_human']}"
            out.append(gg)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "osdr_projection"))
    ap.add_argument("--p", type=float, default=0.01)
    ap.add_argument("--fdr", type=float, default=None)
    ap.add_argument("--top", type=int, default=None)
    ap.add_argument("--orthologs", default=None,
                    help="MGI HOM_MouseHumanSequence.rpt path")
    ap.add_argument("--human", action="store_true",
                    help="genes are already human symbols; no mapping")
    ap.add_argument("--contrast", default=DEFAULT_CONTRAST)
    args = ap.parse_args()

    global OUT
    OUT = Path(args.out)
    OUT.mkdir(exist_ok=True)
    client = SpokeClient()

    genes = load_signature(args.table, args.p, args.fdr, args.top, args.contrast)
    crit = f"FDR < {args.fdr}" if args.fdr is not None else f"p < {args.p}"
    print(f"Spaceflight signature: {len(genes)} mouse genes at {crit}"
          + (f" (top {args.top} by p)" if args.top else "")
          + f" — {sum(1 for g in genes if g['adj_p'] < 0.1)} at FDR < 0.1")
    genes = add_orthologs(genes, args.orthologs, args.human)
    unmapped = sum(1 for g in genes if not g["human"])
    if args.orthologs:
        print(f"MGI ortholog mapping: {len(genes) - unmapped} human targets, "
              f"{unmapped} mouse genes with no human ortholog")
    genes = [g for g in genes if g["human"]] + [g for g in genes if not g["human"]]

    disease_hits = defaultdict(list)   # disease name -> [gene, ...]
    pathway_hits = defaultdict(list)
    disease_meta, pathway_meta = {}, {}
    edges = []

    for i, g in enumerate(genes, 1):
        if not g["human"]:
            g["in_spoke"] = False
            continue
        # exact-name lookup in SPOKE (human gene symbols)
        try:
            found = client.node("Gene", "name", g["human"])
        except Exception as e:
            print(f"  ! {g['human']}: lookup failed ({e})")
            g["in_spoke"] = False
            continue
        g["in_spoke"] = bool(found)
        if not found:
            continue

        try:
            graph = client.neighborhood("Gene", "name", g["human"],
                                        node_filters=["Disease", "Pathway"])
        except Exception as e:
            print(f"  ! {g['human']}: neighborhood failed ({e})")
            continue
        nodes, glinks = split_graph(graph)
        by_id = {n["id"]: n for n in nodes}
        root = next((n["id"] for n in nodes
                     if n["neo4j_type"] == "Gene" and n["properties"].get("name") == g["human"]), None)
        n_dis = n_pw = 0
        for e in glinks:
            other = by_id.get(e["target"] if e["source"] == root else e["source"])
            if other is None:
                continue
            props, name = other["properties"], other["properties"].get("name")
            if other["neo4j_type"] == "Disease" and e["neo4j_type"] == "ASSOCIATES_DaG":
                disease_hits[name].append(g["human"])
                disease_meta[name] = props.get("identifier", "")
                edges.append([g["human"], "ASSOCIATES_DaG", "Disease", name])
                n_dis += 1
            elif other["neo4j_type"] == "Pathway":
                pathway_hits[name].append(g["human"])
                pathway_meta[name] = props.get("identifier", "")
                edges.append([g["human"], e["neo4j_type"], "Pathway", name])
                n_pw += 1
        g["n_diseases"], g["n_pathways"] = n_dis, n_pw
        print(f"  [{i:2d}/{len(genes)}] {g['mouse']:12s} -> {g['human']:12s} "
              f"in SPOKE: {n_dis} diseases, {n_pw} pathways")
        time.sleep(0.2)  # be polite to the public API

    with (OUT / "gene_matches.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mouse_gene", "human_gene", "homology", "log2fc", "p", "adj_p",
                    "in_spoke", "n_diseases", "n_pathways"])
        for g in genes:
            w.writerow([g["mouse"], g["human"], g.get("homology", ""),
                        round(g["log2fc"], 3), g["p"], round(g["adj_p"], 3),
                        g.get("in_spoke", False), g.get("n_diseases", 0), g.get("n_pathways", 0)])

    def dump_ranked(path, hits, meta, kind):
        ranked = sorted(hits.items(), key=lambda kv: -len(set(kv[1])))
        with path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow([kind, "identifier", "n_signature_genes", "genes"])
            for name, gg in ranked:
                gg = sorted(set(gg))
                w.writerow([name, meta.get(name, ""), len(gg), ";".join(gg)])
        return ranked

    dis = dump_ranked(OUT / "diseases.csv", disease_hits, disease_meta, "disease")
    pws = dump_ranked(OUT / "pathways.csv", pathway_hits, pathway_meta, "pathway")

    with (OUT / "edges.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "edge_type", "target_type", "target"])
        w.writerows(edges)

    matched = sum(1 for g in genes if g.get("in_spoke"))
    print(f"\n{matched}/{len(genes)} signature genes found in SPOKE")
    print(f"{len(dis)} diseases touched — top 10 by signature-gene count:")
    for name, gg in dis[:10]:
        gg = sorted(set(gg))
        print(f"  {len(gg):2d}  {name}  ({', '.join(gg[:6])}{'…' if len(gg) > 6 else ''})")
    print(f"{len(pws)} pathways touched — top 5:")
    for name, gg in pws[:5]:
        print(f"  {len(set(gg)):2d}  {name}")


if __name__ == "__main__":
    main()

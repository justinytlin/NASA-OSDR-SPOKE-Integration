#!/usr/bin/env python3
"""Predict candidate drugs and gene-therapy targets from a SPOKE projection.

Logic (hypothesis generation, not clinical evidence):

DRUGS — signature reversal via SPOKE's CMAP/LINCS compound-gene edges:
  a compound scores + when it regulates a signature gene OPPOSITE to its
  spaceflight direction (flight-up gene, compound downregulates it), and
  scores - when concordant (would push the gene further the wrong way).
  score = sum over genes of w_g * (+1 reversal / -1 concordant), w_g = |log2fc|

GENE THERAPY — signature genes ranked as intervention targets:
  priority = |log2fc| * -log10(p) * (1 + n neurodegenerative-disease links).
  Proposed modality: knockdown (RNAi/ASO) for flight-upregulated genes,
  restoration/overexpression (AAV) for flight-downregulated ones. Genes with
  few known compound regulators are flagged as gene-therapy-first targets.

Usage:
  python3 examples/predict_treatments.py PROJECTION_DIR --out OUT_DIR
  (PROJECTION_DIR needs gene_matches.csv + edges.csv from osdr_spoke_projection.py)
"""

import argparse
import csv
import math
import re
import time
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spoke_api import SpokeClient, split_graph

NEURO_RE = re.compile(
    r"parkinson|alzheimer|neurodegen|amyotroph|huntington|dementia|"
    r"central nervous system|nervous system disease|brain disease|"
    r"movement disease|tauopathy|synucleinopathy", re.I)

UP_EDGES = {"UPREGULATES_CuG"}
DOWN_EDGES = {"DOWNREGULATES_CdG"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("projection_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-compounds", type=int, default=20,
                    help="how many top-scoring compounds get a disease-context lookup")
    args = ap.parse_args()
    proj = Path(args.projection_dir)
    out = Path(args.out)
    out.mkdir(exist_ok=True)
    client = SpokeClient()

    genes = [g for g in csv.DictReader(open(proj / "gene_matches.csv"))
             if g["in_spoke"] == "True" and g["human_gene"] not in ("", "NA")]
    for g in genes:
        g["log2fc"] = float(g["log2fc"])
        g["p"] = float(g["p"])
    print(f"{len(genes)} SPOKE-matched signature genes")

    # neurodegenerative-disease links per gene, from the projection's edge list
    neuro_links = defaultdict(list)
    for e in csv.DictReader(open(proj / "edges.csv")):
        if e["target_type"] == "Disease" and NEURO_RE.search(e["target"]):
            neuro_links[e["gene"]].append(e["target"])

    # ── drug scoring: pull compound regulators of each signature gene ──
    comp = {}   # compound name -> record
    for i, g in enumerate(genes, 1):
        gname, flight_up = g["human_gene"], g["log2fc"] > 0
        w = abs(g["log2fc"])
        try:
            graph = client.neighborhood("Gene", "name", gname, node_filters=["Compound"])
        except Exception as exc:
            print(f"  ! {gname}: {exc}")
            continue
        nodes, edges = split_graph(graph)
        by_id = {n["id"]: n for n in nodes}
        root = next((n["id"] for n in nodes
                     if n["neo4j_type"] == "Gene" and n["properties"].get("name") == gname), None)
        for e in edges:
            et = e["neo4j_type"]
            if et not in UP_EDGES | DOWN_EDGES:
                continue
            other = by_id.get(e["target"] if e["source"] == root else e["source"])
            if other is None or other["neo4j_type"] != "Compound":
                continue
            props = other["properties"]
            cname = props.get("name") or props.get("identifier")
            c = comp.setdefault(cname, {
                "name": cname, "identifier": props.get("identifier", ""),
                "max_phase": props.get("max_phase", ""),
                "score": 0.0, "reversed": [], "concordant": []})
            comp_up = et in UP_EDGES
            if comp_up != flight_up:      # compound opposes the flight change
                c["score"] += w
                c["reversed"].append(f"{gname}(cmpd-{'up' if comp_up else 'dn'})")
            else:                          # compound pushes the same way
                c["score"] -= w
                c["concordant"].append(gname)
        if i % 10 == 0:
            print(f"  [{i}/{len(genes)}] scored, {len(comp)} compounds so far")
        time.sleep(0.15)

    ranked = sorted(comp.values(), key=lambda c: -c["score"])
    print(f"{len(comp)} compounds touched; scoring done")

    # ── disease context for the top compounds ──
    for c in ranked[:args.top_compounds]:
        try:
            graph = client.neighborhood("Compound", "name", c["name"],
                                        node_filters=["Disease"])
            nodes, edges = split_graph(graph)
            by_id = {n["id"]: n for n in nodes}
            croot = next((n["id"] for n in nodes if n["neo4j_type"] == "Compound"
                          and n["properties"].get("name") == c["name"]), None)
            treats, neuro = [], []
            for e in edges:
                if e["neo4j_type"] not in ("TREATS_CtD", "IN_CLINICAL_TRIALS_FOR_CictD"):
                    continue
                other = by_id.get(e["target"] if e["source"] == croot else e["source"])
                if other is None:
                    continue
                dname = other["properties"].get("name", "")
                treats.append(dname)
                if NEURO_RE.search(dname):
                    neuro.append(dname)
            c["treats_n"] = len(set(treats))
            c["neuro_indications"] = "; ".join(sorted(set(neuro))[:4])
        except Exception as exc:
            print(f"  ! compound {c['name']}: {exc}")
        time.sleep(0.15)

    with (out / "drug_candidates.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound", "identifier", "max_phase", "reversal_score",
                    "n_genes_reversed", "n_genes_concordant",
                    "known_neuro_indications", "genes_reversed"])
        for c in ranked[:100]:
            w.writerow([c["name"], c["identifier"], c["max_phase"], round(c["score"], 2),
                        len(c["reversed"]), len(c["concordant"]),
                        c.get("neuro_indications", ""), ";".join(c["reversed"][:12])])

    # ── gene-therapy target ranking ──
    rows = []
    for g in genes:
        nd = neuro_links.get(g["human_gene"], [])
        logp = -math.log10(g["p"]) if g["p"] > 0 else 320.0
        priority = abs(g["log2fc"]) * logp * (1 + len(set(nd)))
        n_regulators = sum(1 for c in comp.values()
                           if any(r.startswith(g["human_gene"] + "(") for r in c["reversed"])
                           or g["human_gene"] in c["concordant"])
        rows.append({
            "gene": g["human_gene"], "log2fc": g["log2fc"], "p": g["p"],
            "priority": round(priority, 1),
            "proposed_modality": ("knockdown (RNAi/ASO)" if g["log2fc"] > 0
                                  else "restore/overexpress (AAV)"),
            "n_neuro_diseases": len(set(nd)),
            "neuro_diseases": "; ".join(sorted(set(nd))[:4]),
            "n_compound_regulators": n_regulators,
            "gene_therapy_first": "yes" if n_regulators < 3 else "",
        })
    rows.sort(key=lambda r: -r["priority"])
    with (out / "gene_therapy_targets.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("\nTop 15 drug candidates by reversal score:")
    for c in ranked[:15]:
        print(f"  {c['score']:7.2f}  {str(c['name'])[:38]:38s} phase={c['max_phase'] or '?'} "
              f"rev={len(c['reversed'])} conc={len(c['concordant'])} "
              f"{('NEURO: ' + c['neuro_indications']) if c.get('neuro_indications') else ''}")
    print("\nTop 12 gene-therapy targets:")
    for r in rows[:12]:
        print(f"  {r['priority']:8.1f}  {r['gene']:10s} {r['proposed_modality']:26s} "
              f"neuro-links={r['n_neuro_diseases']} {'[gene-therapy-first]' if r['gene_therapy_first'] else ''}")


if __name__ == "__main__":
    main()

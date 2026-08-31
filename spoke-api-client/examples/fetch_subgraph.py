#!/usr/bin/env python3
"""Fetch a local SPOKE subgraph: 1-hop neighborhoods of a set of genes.

Writes one JSON line per gene: {"gene": ..., "nodes": [...], "edges": [...]}
Usage: fetch_subgraph.py OUT.jsonl GENE_MATCHES.csv [GENE_MATCHES2.csv ...]
"""

import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spoke_api import SpokeClient, split_graph

out_path, match_files = sys.argv[1], sys.argv[2:]
genes = set()
for mf in match_files:
    for g in csv.DictReader(open(mf)):
        if g["in_spoke"] == "True" and g["human_gene"] not in ("", "NA"):
            genes.add(g["human_gene"])
genes = sorted(genes)
print(f"{len(genes)} unique SPOKE-matched genes across {len(match_files)} signatures")

client = SpokeClient()
with open(out_path, "w") as out:
    for i, gname in enumerate(genes, 1):
        try:
            graph = client.neighborhood("Gene", "name", gname)
        except Exception as exc:
            print(f"  ! {gname}: {exc}")
            continue
        nodes, edges = split_graph(graph)
        rec = {
            "gene": gname,
            "nodes": [{"id": n["id"], "type": n["neo4j_type"],
                       "name": n["properties"].get("name")} for n in nodes],
            "edges": [{"s": e["source"], "t": e["target"], "type": e["neo4j_type"]}
                      for e in edges],
        }
        out.write(json.dumps(rec) + "\n")
        if i % 10 == 0:
            print(f"  [{i}/{len(genes)}] fetched")
        time.sleep(0.15)
print("done")

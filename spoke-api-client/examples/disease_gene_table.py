#!/usr/bin/env python3
"""Example workflow: search a disease, pull its neighborhood, and build a
gene-association table — the shape of output you'd join against an external
dataset (e.g. differentially expressed genes from a NASA OSDR study).

Usage:
  python3 examples/disease_gene_table.py "multiple sclerosis"
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spoke_api import SpokeClient, split_graph


def main():
    query = sys.argv[1] if len(sys.argv) > 1 else "multiple sclerosis"
    client = SpokeClient()

    hits = client.search(query, node_type="Disease")
    if not hits:
        sys.exit(f"No Disease nodes match {query!r}")
    top = hits[0]
    print(f"Top hit: {top['name']} ({top['identifier']})")

    graph = client.neighborhood(
        "Disease", "identifier", top["identifier"],
        node_filters=["Gene"],
    )
    nodes, edges = split_graph(graph)
    genes = {n["id"]: n for n in nodes if n.get("neo4j_type") == "Gene"}
    print(f"Neighborhood: {len(nodes)} nodes, {len(edges)} edges, {len(genes)} genes")

    out = Path(f"{top['identifier'].replace(':', '_')}_genes.csv")
    with out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["disease", "disease_id", "edge_type", "gene_symbol", "gene_id"])
        for e in edges:
            gene = genes.get(e["source"]) or genes.get(e["target"])
            if gene:
                props = gene.get("properties", {})
                writer.writerow([
                    top["name"], top["identifier"],
                    e.get("neo4j_type", ""),
                    props.get("name", ""),
                    props.get("identifier", ""),
                ])
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

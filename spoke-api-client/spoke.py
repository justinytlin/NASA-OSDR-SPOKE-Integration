#!/usr/bin/env python3
"""Command-line interface for the SPOKE REST API.

Examples:
  ./spoke.py version
  ./spoke.py types
  ./spoke.py search "multiple sclerosis" --type Disease
  ./spoke.py node Disease identifier DOID:2377
  ./spoke.py neighborhood Disease identifier DOID:2377 --output ms.json
  ./spoke.py neighborhood Compound identifier DB01611 --nodes Gene Protein --format csv
"""

import argparse
import csv
import io
import json
import sys

from spoke_api import SpokeClient, split_graph


def emit(data, args, default_name: str):
    """Write result as JSON or CSV to --output or stdout."""
    if args.format == "csv":
        rows = data
        if isinstance(data, dict):  # graph response -> node/edge tables
            nodes, edges = split_graph(data)
            rows = nodes + edges
        if not isinstance(rows, list) or not rows:
            sys.exit("No tabular rows to write as CSV.")
        rows = [r.get("data", r) if isinstance(r, dict) else {"value": r} for r in rows]
        # flatten nested "properties" dicts (graph responses) into columns
        rows = [
            {**{k: v for k, v in r.items() if k != "properties"}, **r.get("properties", {})}
            if isinstance(r.get("properties"), dict) else r
            for r in rows
        ]
        fields: list[str] = []
        for r in rows:
            fields += [k for k in r if k not in fields]
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        text = buf.getvalue()
    else:
        text = json.dumps(data, indent=2)

    if args.output:
        path = args.output
        if "." not in path.rsplit("/", 1)[-1]:
            path += f".{args.format}"
        with open(path, "w") as f:
            f.write(text)
        print(f"Wrote {path}", file=sys.stderr)
    else:
        print(text)


def main():
    parser = argparse.ArgumentParser(description="Query the SPOKE knowledge graph REST API")
    parser.add_argument("--format", choices=["json", "csv"], default="json")
    parser.add_argument("--output", "-o", help="write to file instead of stdout")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("version", help="SPOKE version and data-source update times")
    sub.add_parser("metagraph", help="node/edge type metagraph")
    sub.add_parser("types", help="node types, edge types, default query params")

    p = sub.add_parser("search", help="search nodes by name (lucene syntax)")
    p.add_argument("query")
    p.add_argument("--type", dest="node_type", help="restrict to a node type, e.g. Disease")

    p = sub.add_parser("node", help="fetch one node: node <type> <attribute> <value>")
    p.add_argument("node_type")
    p.add_argument("attribute")
    p.add_argument("value")

    p = sub.add_parser("neighborhood", help="neighborhood graph around a node")
    p.add_argument("node_type")
    p.add_argument("attribute")
    p.add_argument("value")
    p.add_argument("--depth", type=int)
    p.add_argument("--nodes", nargs="*", dest="node_filters", metavar="TYPE",
                   help="only include these node types")
    p.add_argument("--edges", nargs="*", dest="edge_filters", metavar="TYPE",
                   help="only include these edge types")
    p.add_argument("--cutoff", nargs=2, action="append", default=[],
                   metavar=("NAME", "VALUE"),
                   help="e.g. --cutoff cutoff_DaG_textmining 3.0")

    p = sub.add_parser("expand", help="expand a node by internal id")
    p.add_argument("node_type")
    p.add_argument("node_id", type=int)

    p = sub.add_parser("sea", help="SEA compound similarity by SMILES or ZINC id")
    p.add_argument("smiles_or_zinc")

    args = parser.parse_args()
    client = SpokeClient()

    if args.command == "version":
        result = client.version()
    elif args.command == "metagraph":
        result = client.metagraph()
    elif args.command == "types":
        result = client.types()
    elif args.command == "search":
        result = client.search(args.query, args.node_type)
    elif args.command == "node":
        result = client.node(args.node_type, args.attribute, args.value)
    elif args.command == "neighborhood":
        cutoffs = {name: value for name, value in args.cutoff}
        result = client.neighborhood(
            args.node_type, args.attribute, args.value,
            depth=args.depth, node_filters=args.node_filters,
            edge_filters=args.edge_filters, **cutoffs,
        )
    elif args.command == "expand":
        result = client.expand(args.node_type, args.node_id)
    elif args.command == "sea":
        result = client.sea(args.smiles_or_zinc)

    emit(result, args, args.command)


if __name__ == "__main__":
    main()

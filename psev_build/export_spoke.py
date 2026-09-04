#!/usr/bin/env python3
"""Export the node table and edge list of a SPOKE Neo4j instance for PSEV building.

Reads Neo4j credentials from environment variables or a .env file (searched in
the current directory, this script's directory, and its parent):

    NEO4J_URI=bolt://localhost:7687      (default)
    NEO4J_USER=neo4j
    NEO4J_PASSWORD=...
    NEO4J_DATABASE=neo4j                 (default; SPOKE deployments often use "spoke")

(The spoke-cli names KNOWLEDGE_GRAPH_URI / _USERNAME / _PASSWORD / _DATABASE are
accepted as aliases.)

Usage:
  export_spoke.py stats  [--out DIR]
      Count nodes per label and relationships per type -> DIR/graph_stats.json
  export_spoke.py nodes  [--out DIR] [--include-types A,B | --exclude-types A,B]
                         [--protein all|human|human-reviewed]
      -> DIR/node_info.tsv  (Node, Node_Name, Node_Type, Node_Index) in the same
         shape as SPOKE_NASA/spoke_v_2/node_info.tsv, plus DIR/node_internal_ids.npy
         (Neo4j internal id per Node_Index, used to join the edge export)
  export_spoke.py edges  [--out DIR] [--exclude-rel-types X,Y] [--rel-types X,Y]
      -> DIR/edges/<REL_TYPE>.npy  (int64 array, shape (n, 2), internal ids of
         (start, end)); one file per relationship type, skipped if it already
         exists so an interrupted export can be resumed.

Node identity follows the v2 convention: Node = n.identifier (Entrez id for
Gene, UniProt for Protein, DOID for Disease, ...), Node_Name = n.name,
Node_Type = first label.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def load_dotenv():
    for d in (Path.cwd(), HERE, HERE.parent):
        p = d / ".env"
        if p.exists():
            for line in p.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))
            return p
    return None


def get_driver():
    from neo4j import GraphDatabase
    env = load_dotenv()
    uri = os.environ.get("NEO4J_URI") or os.environ.get("KNOWLEDGE_GRAPH_URI") or "bolt://localhost:7687"
    user = os.environ.get("NEO4J_USER") or os.environ.get("KNOWLEDGE_GRAPH_USERNAME") or "neo4j"
    pw = os.environ.get("NEO4J_PASSWORD") or os.environ.get("KNOWLEDGE_GRAPH_PASSWORD")
    db = os.environ.get("NEO4J_DATABASE") or os.environ.get("KNOWLEDGE_GRAPH_DATABASE") or "neo4j"
    if not pw:
        sys.exit("No Neo4j password found. Put NEO4J_USER / NEO4J_PASSWORD "
                 "(optionally NEO4J_URI, NEO4J_DATABASE) in a .env file next to "
                 "this script or export them as environment variables.")
    print(f"connecting to {uri} db={db} as {user}" + (f" (.env: {env})" if env else ""))
    driver = GraphDatabase.driver(uri, auth=(user, pw), connection_timeout=30,
                                  max_connection_lifetime=3600 * 6)
    driver.verify_connectivity()
    # SPOKE servers usually name the database "spoke"; fall back to the server's
    # default database if the configured one does not exist.
    with driver.session(database="system") as s:
        dbs = [(r["name"], r.get("default"), r.get("type")) for r in s.run("SHOW DATABASES")]
    names = [d[0] for d in dbs if d[2] != "system"]
    if db not in names:
        default = next((d[0] for d in dbs if d[1]), None) or (names[0] if names else db)
        print(f"database '{db}' not found; available {names}; using '{default}'")
        db = default
    return driver, db


def internal_id_expr(var):
    # elementId() is "<n>:<db-uuid>:<internal id>" in Neo4j 5+/2025+; the last
    # segment is the numeric internal id. Falls back to id() on old servers.
    return f"toInteger(split(elementId({var}), ':')[-1])"


def run_stats(driver, db, out):
    stats = {"labels": {}, "relationship_types": {}}
    with driver.session(database=db) as s:
        labels = [r["label"] for r in s.run("CALL db.labels() YIELD label RETURN label")]
        for lab in labels:
            n = s.run(f"MATCH (n:`{lab}`) RETURN count(n) AS c").single()["c"]
            stats["labels"][lab] = n
            print(f"  {lab:24s} {n:>12,}")
        rels = [r["relationshipType"] for r in
                s.run("CALL db.relationshipTypes() YIELD relationshipType RETURN relationshipType")]
        for rt in rels:
            n = s.run(f"MATCH ()-[r:`{rt}`]->() RETURN count(r) AS c").single()["c"]
            stats["relationship_types"][rt] = n
            print(f"  {rt:32s} {n:>12,}")
    stats["total_nodes"] = sum(stats["labels"].values())
    stats["total_edges"] = sum(stats["relationship_types"].values())
    stats["neo4j_version"] = driver.get_server_info().agent
    stats["exported"] = time.strftime("%Y-%m-%d %H:%M:%S")
    out.mkdir(parents=True, exist_ok=True)
    (out / "graph_stats.json").write_text(json.dumps(stats, indent=2))
    print(f"total nodes {stats['total_nodes']:,}  total edges {stats['total_edges']:,}")
    print(f"wrote {out / 'graph_stats.json'}")


def run_nodes(driver, db, out, include, exclude, protein="all"):
    out.mkdir(parents=True, exist_ok=True)
    stats_file = out / "graph_stats.json"
    with driver.session(database=db) as s:
        labels = [r["label"] for r in s.run("CALL db.labels() YIELD label RETURN label")]
        if include:
            missing = set(include) - set(labels)
            if missing:
                sys.exit(f"unknown labels: {sorted(missing)}; available: {labels}")
            labels = [l for l in labels if l in include]
        if exclude:
            labels = [l for l in labels if l not in exclude]
        labels = sorted(labels)
        print("exporting labels:", labels)
        frames = []
        for lab in labels:
            t0 = time.time()
            if lab == "Protein" and protein != "all":
                # human proteins = encoded by the Homo sapiens Organism node
                # (Protein nodes carry org_ncbi_id/org_name, but the edge is indexed)
                # p.reviewed is a string: 'Reviewed, From SwissProt' / 'Unreviewed, From TrEMBL'
                cond = " AND n.reviewed STARTS WITH 'Reviewed'" if protein == "human-reviewed" else ""
                q = (f"MATCH (:Organism {{identifier:'9606'}})-[:ENCODES_OeP]->(n:Protein) "
                     f"WHERE true{cond} RETURN DISTINCT {internal_id_expr('n')} AS iid, "
                     f"n.identifier AS identifier, n.name AS name")
                print(f"  Protein filter: {protein}")
            else:
                q = (f"MATCH (n:`{lab}`) RETURN {internal_id_expr('n')} AS iid, "
                     f"n.identifier AS identifier, n.name AS name")
            rows = []
            for r in s.run(q):
                rows.append((r["iid"], r["identifier"], r["name"]))
            df = pd.DataFrame(rows, columns=["iid", "Node", "Node_Name"])
            df["Node_Type"] = lab
            frames.append(df)
            print(f"  {lab:24s} {len(df):>12,}  ({time.time() - t0:.0f}s)")
    nodes = pd.concat(frames, ignore_index=True)
    # sort like v2 did (by type, then identifier as string) for stable indices
    nodes["Node"] = nodes["Node"].astype(str)
    nodes["Node_Name"] = nodes["Node_Name"].fillna(nodes["Node"]).astype(str)
    nodes = nodes.sort_values(["Node_Type", "Node"]).reset_index(drop=True)
    dup = nodes.duplicated(["Node", "Node_Type"]).sum()
    if dup:
        print(f"WARNING: {dup} duplicate (identifier, type) pairs; keeping first")
        nodes = nodes.drop_duplicates(["Node", "Node_Type"]).reset_index(drop=True)
    nodes["Node_Index"] = np.arange(len(nodes))
    nodes[["Node", "Node_Name", "Node_Type", "Node_Index"]].to_csv(
        out / "node_info.tsv", sep="\t", index=False)
    np.save(out / "node_internal_ids.npy", nodes["iid"].values.astype(np.int64))
    print(f"wrote {out / 'node_info.tsv'} ({len(nodes):,} nodes)")
    print(nodes.Node_Type.value_counts().to_string())


def run_edges(driver, db, out, rel_types, exclude):
    edir = out / "edges"
    edir.mkdir(parents=True, exist_ok=True)
    with driver.session(database=db) as s:
        rels = [r["relationshipType"] for r in
                s.run("CALL db.relationshipTypes() YIELD relationshipType RETURN relationshipType")]
        if rel_types:
            rels = [r for r in rels if r in rel_types]
        if exclude:
            rels = [r for r in rels if r not in exclude]
        rels = sorted(rels)
        print(f"exporting {len(rels)} relationship types")
        for rt in rels:
            target = edir / f"{rt}.npy"
            if target.exists():
                print(f"  {rt:32s} exists, skipping")
                continue
            t0 = time.time()
            q = (f"MATCH (a)-[r:`{rt}`]->(b) RETURN {internal_id_expr('a')} AS s, "
                 f"{internal_id_expr('b')} AS t")
            chunks, buf, n = [], [], 0
            for r in s.run(q):
                buf.append((r["s"], r["t"]))
                if len(buf) >= 1_000_000:
                    chunks.append(np.array(buf, dtype=np.int64))
                    n += len(buf)
                    buf = []
                    print(f"    {rt}: {n:,} rows ({time.time() - t0:.0f}s)", flush=True)
            if buf:
                chunks.append(np.array(buf, dtype=np.int64))
                n += len(buf)
            arr = np.concatenate(chunks) if chunks else np.zeros((0, 2), dtype=np.int64)
            tmp = target.with_suffix(".tmp.npy")
            np.save(tmp, arr)
            tmp.rename(target)
            print(f"  {rt:32s} {n:>12,}  ({time.time() - t0:.0f}s)", flush=True)
    print(f"edge files in {edir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["stats", "nodes", "edges"])
    ap.add_argument("--out", default=str(HERE / "spoke_current"))
    ap.add_argument("--include-types", default=None, help="comma-separated node labels to keep")
    ap.add_argument("--exclude-types", default=None, help="comma-separated node labels to drop")
    ap.add_argument("--protein", default="all", choices=["all", "human", "human-reviewed"],
                    help="nodes: keep all Protein nodes, only human (Organism 9606), "
                         "or only human Swiss-Prot (reviewed) proteins")
    ap.add_argument("--rel-types", default=None, help="comma-separated relationship types to export")
    ap.add_argument("--exclude-rel-types", default=None, help="comma-separated relationship types to skip")
    a = ap.parse_args()
    split = lambda s: set(x.strip() for x in s.split(",") if x.strip()) if s else None
    driver, db = get_driver()
    out = Path(a.out)
    try:
        if a.cmd == "stats":
            run_stats(driver, db, out)
        elif a.cmd == "nodes":
            run_nodes(driver, db, out, split(a.include_types), split(a.exclude_types), a.protein)
        else:
            run_edges(driver, db, out, split(a.rel_types), split(a.exclude_rel_types))
    finally:
        driver.close()


if __name__ == "__main__":
    main()

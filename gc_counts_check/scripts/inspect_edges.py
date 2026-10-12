#!/usr/bin/env python3
"""Inspect Anatomy–Gene edge properties in the live SPOKE Neo4j (needs the SSH tunnel and psev_build/.env).

Reports, for EXPRESSES_AeG, UPREGULATES_AuG and DOWNREGULATES_AdG: the property keys, five sample
edges with their properties, the sources, and the number of distinct Anatomy nodes and genes per type.
Nothing is written to the database. Usage: inspect_edges.py
"""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "psev_build"))
import export_spoke  # noqa: E402  (reads .env, provides the driver)


def main():
    export_spoke.load_dotenv()
    drv, db = export_spoke.get_driver()
    with drv.session(database=db) as s:
        for rel in ["EXPRESSES_AeG", "UPREGULATES_AuG", "DOWNREGULATES_AdG"]:
            print(f"\n=== {rel}")
            keys = s.run(f"MATCH ()-[r:{rel}]->() WITH r LIMIT 2000 UNWIND keys(r) AS k RETURN k, count(*) AS n ORDER BY n DESC").data()
            print("property keys (first 2000 edges):", {d['k']: d['n'] for d in keys})
            for d in s.run(f"MATCH (a)-[r:{rel}]->(g) RETURN a.name AS anatomy, g.name AS gene, properties(r) AS p LIMIT 5").data():
                print("  ", d["anatomy"], "->", d["gene"], d["p"])
            src = s.run(f"MATCH ()-[r:{rel}]->() WITH r LIMIT 20000 RETURN r.sources AS src, count(*) AS n ORDER BY n DESC LIMIT 8").data()
            print("sources (first 20k edges):", [(d["src"], d["n"]) for d in src])
            cnt = s.run(f"MATCH (a:Anatomy)-[r:{rel}]->(g:Gene) RETURN count(DISTINCT a) AS anatomies, count(DISTINCT g) AS genes, count(r) AS edges").data()[0]
            print("coverage:", cnt)
            per = s.run(f"MATCH (a:Anatomy)-[r:{rel}]->(:Gene) RETURN a.name AS a, count(r) AS n ORDER BY n DESC LIMIT 6").data()
            print("most-linked anatomy nodes:", [(d["a"], d["n"]) for d in per])
    drv.close()


if __name__ == "__main__":
    main()

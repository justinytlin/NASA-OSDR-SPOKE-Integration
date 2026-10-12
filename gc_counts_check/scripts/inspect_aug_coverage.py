#!/usr/bin/env python3
"""Per-Anatomy counts of UPREGULATES_AuG / DOWNREGULATES_AdG / EXPRESSES_AeG edges and the AuG log2fc
distribution, from the live SPOKE (read-only). Writes results/anatomy_gene_edge_counts.tsv."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "psev_build"))
import export_spoke
export_spoke.load_dotenv(); drv, db = export_spoke.get_driver()
with drv.session(database=db) as s:
    q = """MATCH (a:Anatomy) WHERE (a)-[:UPREGULATES_AuG|DOWNREGULATES_AdG|EXPRESSES_AeG]->(:Gene)
           RETURN a.identifier AS id, a.name AS name,
             size([(a)-[:UPREGULATES_AuG]->(:Gene) | 1]) AS n_up,
             size([(a)-[r:UPREGULATES_AuG]->(:Gene) WHERE r.log2fc >= 2 | 1]) AS n_up_log2fc_ge2,
             size([(a)-[:DOWNREGULATES_AdG]->(:Gene) | 1]) AS n_down,
             size([(a)-[:EXPRESSES_AeG]->(:Gene) | 1]) AS n_expr
           ORDER BY n_up DESC"""
    df = pd.DataFrame(s.run(q).data())
    dist = s.run("MATCH ()-[r:UPREGULATES_AuG]->() RETURN percentileCont(r.log2fc, 0.1) AS p10, percentileCont(r.log2fc, 0.5) AS p50, "
                 "percentileCont(r.log2fc, 0.9) AS p90, min(r.log2fc) AS mn, max(r.log2fc) AS mx, count(r) AS n").data()[0]
drv.close()
out = Path(__file__).resolve().parent.parent / "results/anatomy_gene_edge_counts.tsv"; df.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
print("AuG log2fc distribution:", {k: round(v, 2) if isinstance(v, float) else v for k, v in dist.items()})
print(f"{len(df)} Anatomy nodes with any Anatomy-Gene edge; with >=1 AuG: {(df.n_up > 0).sum()}, >=50 AuG: {(df.n_up >= 50).sum()}, >=200: {(df.n_up >= 200).sum()}")
print(df.to_string(index=False))

"""Mouse -> human ortholog mapping from the MGI/JAX homology report.

Source: https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
Rows sharing a "DB Class Key" form one homology class containing the mouse
and human members of that ortholog group.
"""

import csv
from collections import defaultdict


def load_mgi(path):
    """Parse the MGI report into a mouse->human ortholog map.

    Returns (by_entrez, by_symbol): each maps a mouse key to a list of
    {"symbol", "entrez", "n_mouse", "n_human"} human-gene dicts, where the
    counts describe the homology class (1/1 = clean 1:1 ortholog).
    """
    classes = defaultdict(lambda: {"mouse": [], "human": []})
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = row["DB Class Key"]
            org = row["Common Organism Name"]
            entry = {"symbol": row["Symbol"], "entrez": row["EntrezGene ID"]}
            if org.startswith("mouse"):
                classes[key]["mouse"].append(entry)
            elif org == "human":
                classes[key]["human"].append(entry)

    by_entrez, by_symbol = {}, {}
    for c in classes.values():
        if not c["mouse"] or not c["human"]:
            continue
        humans = [{"symbol": h["symbol"], "entrez": h["entrez"],
                   "n_mouse": len(c["mouse"]), "n_human": len(c["human"])}
                  for h in c["human"]]
        for m in c["mouse"]:
            if m["entrez"]:
                by_entrez[m["entrez"]] = humans
            by_symbol[m["symbol"]] = humans
    return by_entrez, by_symbol


def map_gene(mouse_entrez, mouse_symbol, by_entrez, by_symbol):
    """Human orthologs for one mouse gene; entrez match wins over symbol."""
    if mouse_entrez and str(mouse_entrez) in by_entrez:
        return by_entrez[str(mouse_entrez)]
    return by_symbol.get(mouse_symbol, [])

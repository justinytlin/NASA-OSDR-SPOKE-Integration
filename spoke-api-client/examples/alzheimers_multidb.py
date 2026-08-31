#!/usr/bin/env python3
"""Multi-database view of Alzheimer's disease (DOID:10652) from SPOKE.

Pulls the disease neighborhood and splits it into per-source tables:
  genes         ASSOCIATES_DaG   (GWAS Catalog, DISEASES, OMIM)
  drugs         TREATS_CtD       (ChEMBL / DrugCentral, with clinical phase)
  trials        IN_CLINICAL_TRIALS_FOR_CictD (ClinicalTrials.gov via ChEMBL)
  contra        CONTRAINDICATES_CcD (DrugCentral)
  proteins      INCREASEDIN_PiD  (proteomics / Cell Taxonomy)
  prevalence    PREVALENCE_DpL   (IHME Global Burden of Disease, by location)
Then a 2-hop hop through APOE for non-disease databases:
  pathways      PARTICIPATES_GpPW (KEGG / Reactome / WikiPathways)
  go terms      PARTICIPATES_GpBP/GpMF/GpCC (Gene Ontology)
  anatomy       EXPRESSES/UPREGULATES/DOWNREGULATES (BGee)

Writes one CSV per table into ./alzheimers_out/.
"""

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spoke_api import SpokeClient, split_graph

AD_ID = "DOID:10652"
OUT = Path(__file__).resolve().parents[1] / "alzheimers_out"


def write_csv(name: str, header: list[str], rows: list[list]):
    path = OUT / f"{name}.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path.name}: {len(rows)} rows")


def index_graph(graph):
    nodes, edges = split_graph(graph)
    return {n["id"]: n for n in nodes}, edges


def neighbor_of(root_id, edge, by_id):
    """The endpoint of the edge that is not the query's root node."""
    other = edge["target"] if edge["source"] == root_id else edge["source"]
    return by_id.get(other)


def main():
    OUT.mkdir(exist_ok=True)
    client = SpokeClient()

    ad = client.node("Disease", "identifier", AD_ID)[0]
    print(f"Disease: {ad['name']} ({AD_ID}), sources: {ad.get('source')}")

    # ── 1-hop neighborhood, default cutoffs ──────────────────────────────
    graph = client.neighborhood("Disease", "identifier", AD_ID)
    by_id, edges = index_graph(graph)
    root_id = next(i for i, n in by_id.items()
                   if n["properties"].get("identifier") == AD_ID)
    print(f"1-hop neighborhood: {len(by_id)} nodes, {len(edges)} edges")

    def rows_for(edge_type, cols):
        """One row per edge of edge_type: neighbor properties + edge properties."""
        out = []
        for e in edges:
            if e["neo4j_type"] != edge_type:
                continue
            nb = neighbor_of(root_id, e, by_id)
            nbp, ep = nb["properties"], e.get("properties", {})
            out.append([nbp.get(c) if src == "n" else ep.get(c) for src, c in cols])
        return out

    print("\nPer-database tables:")
    write_csv("genes_associates", ["gene", "entrez_id", "ensembl", "sources", "gwas_pvalue"],
              rows_for("ASSOCIATES_DaG", [("n", "name"), ("n", "identifier"),
                                          ("n", "ensembl"), ("e", "sources"), ("e", "gwas_pvalue")]))
    write_csv("drugs_treats", ["compound", "chembl_id", "phase", "sources"],
              rows_for("TREATS_CtD", [("n", "name"), ("n", "identifier"),
                                      ("e", "phase"), ("e", "sources")]))
    write_csv("drugs_in_trials", ["compound", "chembl_id", "max_phase"],
              rows_for("IN_CLINICAL_TRIALS_FOR_CictD",
                       [("n", "name"), ("n", "identifier"), ("n", "max_phase")]))
    write_csv("proteins_increased", ["protein", "uniprot", "sources"],
              rows_for("INCREASEDIN_PiD", [("n", "name"), ("n", "identifier"), ("e", "sources")]))
    write_csv("prevalence_by_location",
              ["location", "geonames_id", "percent", "lower", "upper", "year"],
              rows_for("PREVALENCE_DpL", [("n", "name"), ("n", "identifier"),
                                          ("e", "value"), ("e", "lower"),
                                          ("e", "upper"), ("e", "year")]))

    # contraindications need an explicit edge filter
    cgraph = client.neighborhood("Disease", "identifier", AD_ID,
                                 edge_filters=["CONTRAINDICATES_CcD"])
    cby, cedges = index_graph(cgraph)
    croot = next(i for i, n in cby.items()
                 if n["properties"].get("identifier") == AD_ID)
    rows = []
    for e in cedges:
        nb = neighbor_of(croot, e, cby)
        rows.append([nb["properties"].get("name"), nb["properties"].get("identifier")])
    write_csv("drugs_contraindicated", ["compound", "chembl_id"], rows)

    # ── 2-hop: APOE's non-disease context ────────────────────────────────
    apoe = client.search("APOE", node_type="Gene")[0]
    print(f"\n2-hop via gene: {apoe['name']}")
    # Gene identifiers are stored as integers; matching by name is reliable
    ggraph = client.neighborhood("Gene", "name", apoe["name"])
    gby, gedges = index_graph(ggraph)
    groot = next(i for i, n in gby.items()
                 if n["neo4j_type"] == "Gene" and n["properties"].get("name") == apoe["name"])
    print(f"APOE neighborhood: {len(gby)} nodes, {len(gedges)} edges")
    print("  edge types:", dict(Counter(e["neo4j_type"] for e in gedges)))

    def gene_rows(edge_types):
        out = []
        for e in gedges:
            if e["neo4j_type"] in edge_types:
                nb = neighbor_of(groot, e, gby)
                nbp = nb["properties"]
                out.append([e["neo4j_type"], nb["neo4j_type"],
                            nbp.get("name"), nbp.get("identifier"), nbp.get("source")])
        return out

    write_csv("apoe_pathways", ["edge", "type", "name", "identifier", "source"],
              gene_rows({"PARTICIPATES_GpPW"}))
    write_csv("apoe_go_terms", ["edge", "type", "name", "identifier", "source"],
              gene_rows({"PARTICIPATES_GpBP", "PARTICIPATES_GpMF", "PARTICIPATES_GpCC"}))
    write_csv("apoe_celltype_expression", ["edge", "type", "name", "identifier", "source"],
              gene_rows({"EXPRESSEDIN_GeiCT"}))
    write_csv("apoe_proteins_variants", ["edge", "type", "name", "identifier", "source"],
              gene_rows({"ENCODES_GeP", "ASSOCIATES_GaS"}))


if __name__ == "__main__":
    main()

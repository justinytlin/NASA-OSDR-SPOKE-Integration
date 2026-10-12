#!/usr/bin/env python3
"""Cell-type route: rank the CellType nodes by baseline-corrected score for every dataset of a build run.

For tissues with no gene-linked Anatomy node (thymus, retina, bone marrow, placenta, ...) the 2025 graph
still has CellType nodes with Gene edges (59 of them, e.g. retinal rod cell, T cell, erythroid progenitor
cell, syncytiotrophoblast cell). This ranks all CellType nodes of the build per sample (corrected score,
1 = top), pools over samples, and reports the best-ranked cell type whose name matches the tissue, plus
the top five cell types overall. Writes results/celltype_check_<tag>.tsv.
Usage: celltype_check.py --tag 2025upAuG --build /Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only
"""
import argparse, json, re
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent.parent
FAM = {"thymus": r"thymocyte|thymic|\bT cell|NK T cell", "retina": r"photorecept|retina|rod cell|cone cell|horizontal cell|amacrine|bipolar",
       "bone marrow": r"hematopoietic|haematopoietic|erythroid|erythroblast|myeloid|myeloblast|megakaryocyte|granulocyte|neutrophil|monocyte",
       "heart": r"cardiac|cardiomyocyte", "kidney": r"kidney|renal|podocyte|nephron|collecting duct|tubul",
       "colon": r"enterocyte|colon|intestin|goblet|Paneth", "cecum": r"enterocyte|colon|intestin|goblet|Paneth",
       "mammary gland": r"mammary|luminal|lactocyte|breast", "placenta": r"trophoblast|decidu|Hofbauer|placent",
       "eye": r"lens|cornea|retina|photorecept|rod cell|cone cell", "bone": r"osteo|chondro", "mandible": r"osteo|chondro|odonto|amelo",
       "temporal bone": r"osteo|chondro", "optic nerve": r"oligodendro|astrocyte|Schwann|retinal ganglion",
       "liver": r"hepatocyte|hepatic|Kupffer|cholangio", "spleen": r"\bB cell|splen|erythro|macrophage|plasma cell", "lung": r"pneumocyte|alveolar|club cell|lung|bronch|ciliated",
       "skin": r"keratinocyte|melanocyte|epiderm|fibroblast of skin", "hippocampus": r"neuron|pyramidal|granule|astrocyte|oligodendro",
       "cerebellum": r"Purkinje|granule|cerebell|astrocyte|oligodendro", "cerebrum": r"neuron|pyramidal|astrocyte|oligodendro|interneuron",
       "brain": r"neuron|astrocyte|oligodendro|microglia", "adrenal gland": r"adrenal|chromaffin|cortical cell|steroid",
       "adipose tissue": r"adipocyte|preadipocyte|fat cell", "whole blood": r"erythro|neutrophil|monocyte|lymphocyte|\bT cell|\bB cell|NK",
       "skeletal muscle": r"myocyte|myotube|skeletal muscle|satellite|myoblast"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--build", required=True); a = ap.parse_args()
    ni = pd.read_csv(Path(a.build) / "node_info.tsv", sep="\t"); ct = (ni.Node_Type == "CellType").values
    ctn = ni.Node_Name.astype(str).values[ct]
    rows = []
    for e in json.load(open(ROOT / "data/_catalog" / f"manifest_{a.tag}.json")):
        f = ROOT / "psev_out" / f"{e['name']}_GC_ranks.npz"
        if not f.exists(): continue
        z = np.load(f)
        corr = z["scores"].astype(np.float64) - z["mean_weight"][:, None] * z["uniform"].astype(np.float64)[None, :]
        c = corr[:, ct]; r = (-c).argsort(axis=1).argsort(axis=1) + 1; pooled = r.mean(axis=0); pr = pooled.argsort().argsort() + 1
        rx = FAM.get(e["family"]) or FAM.get(e["tissue"])
        m = np.array([bool(re.search(rx, n, re.I)) for n in ctn]) if rx else np.zeros(len(ctn), bool)
        rows.append(dict(dataset=e["base"], tissue=e["tissue"], control_type=e["control_type"], n_samples=len(z["samples"]),
                         best_family_celltype=ctn[m][pr[m].argmin()] if m.any() else None, celltype_rank=int(pr[m].min()) if m.any() else None,
                         n_celltypes=int(ct.sum()), frac_samples_family_top5=float((r[:, m].min(axis=1) <= 5).mean()) if m.any() else None,
                         top5=", ".join(ctn[np.argsort(pooled)[:5]])))
    d = pd.DataFrame(rows); out = ROOT / "results" / f"celltype_check_{a.tag}.tsv"; d.to_csv(out, sep="\t", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 80)
    print(d.sort_values(["tissue", "celltype_rank"]).to_string(index=False)); print("wrote", out)


if __name__ == "__main__":
    main()

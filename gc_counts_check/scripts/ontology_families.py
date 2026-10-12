#!/usr/bin/env python3
"""Tissue families from the UBERON ontology instead of name matching.

A family is a root term plus everything below it through is_a and part_of
(transitively), intersected with the Anatomy nodes of SPOKE v2. Developmental
relations (develops_from) are not followed, so primordia and precursors are
excluded unless UBERON files them as parts of the organ.

ROOTS maps the family keys used by check_anatomy.py to UBERON root term names;
names are resolved against the ontology at load time and a missing name is an error.
"""
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OBO = ROOT / "data/_ontology/uberon-basic.obo"

ROOTS = {
    "hippocampus": ["hippocampal formation"], "thymus": ["thymus"], "spleen": ["spleen"], "liver": ["liver"],
    "liver parenchyma": ["liver parenchyma"], "cerebral cortex": ["cerebral cortex"], "heart": ["heart"], "kidney": ["kidney"],
    "lung": ["lung"], "skeletal muscle": ["skeletal muscle tissue", "skeletal muscle organ, vertebrate"],
    "bone marrow": ["bone marrow"], "skin": ["skin of body"], "intestine": ["intestine"], "testis": ["testis"],
    "blood": ["blood"], "adrenal gland": ["adrenal gland"], "eye": ["camera-type eye"], "retina": ["retina"],
    "brain": ["brain"], "cerebellum": ["cerebellum"], "cerebrum": ["telencephalon"], "optic nerve": ["cranial nerve II"],
    "colon": ["colon"], "cecum": ["caecum"], "mammary gland": ["mammary gland"], "temporal bone": ["temporal bone"],
    "mandible": ["mandible"], "bone tissue": ["bone tissue"], "adipose tissue": ["adipose tissue"],
    "whole blood": ["blood"], "placenta": ["placenta"],
}


def load_obo(path=OBO):
    name, parents, obsolete = {}, defaultdict(set), set()
    tid = None
    for line in open(path, encoding="utf8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            tid = "?"
        elif line.startswith("[") and line.endswith("]"):
            tid = None
        elif tid is not None:
            if line.startswith("id: "):
                tid = line[4:]
            elif line.startswith("name: "):
                name[tid] = line[6:]
            elif line.startswith("is_a: "):
                parents[tid].add(line[6:].split(" ! ")[0].strip())
            elif line.startswith("relationship: part_of "):
                parents[tid].add(line.split()[2])
            elif line.startswith("is_obsolete: true"):
                obsolete.add(tid)
    return name, parents, obsolete


def build_families(spoke_ids):
    """spoke_ids: iterable of SPOKE Anatomy node ids (UBERON:xxxxxxx). Returns {family: set(ids)} and {family: root ids}."""
    name, parents, obsolete = load_obo()
    children = defaultdict(set)
    for c, ps in parents.items():
        for p in ps:
            children[p].add(c)
    by_name = {}
    for t, n in name.items():
        if t.startswith("UBERON:") and t not in obsolete:
            by_name.setdefault(n.lower(), t)
    spoke = set(spoke_ids)
    fams, roots = {}, {}
    for fam, rnames in ROOTS.items():
        rids = []
        for rn in rnames:
            if rn.lower() not in by_name:
                raise KeyError(f"UBERON has no term named '{rn}' (family {fam})")
            rids.append(by_name[rn.lower()])
        seen, stack = set(rids), list(rids)
        while stack:
            for c in children[stack.pop()]:
                if c not in seen:
                    seen.add(c); stack.append(c)
        fams[fam] = seen & spoke; roots[fam] = rids
    return fams, roots


if __name__ == "__main__":
    import pandas as pd, sys, warnings
    warnings.filterwarnings("ignore")
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from families import family_masks
    ni = pd.read_csv(ROOT.parent / "SPOKE_NASA/spoke_v_2/node_info.tsv", sep="\t"); an = ni[ni.Node_Type == "Anatomy"].reset_index(drop=True)
    fams, roots = build_families(an.Node.values)
    nm = family_masks(an.Node_Name)
    id2name = dict(zip(an.Node, an.Node_Name.astype(str)))
    print(f"{'family':18s} {'roots':28s} {'ontology':>8s} {'by name':>8s} {'both':>6s}  examples only in the name-based family")
    for f in ROOTS:
        byname = set(an.Node[nm[f]]) if f in nm else set()
        only_name = sorted(id2name[i] for i in byname - fams[f])[:4]
        print(f"{f:18s} {','.join(roots[f])[:28]:28s} {len(fams[f]):8d} {len(byname):8d} {len(byname & fams[f]):6d}  {only_name}")

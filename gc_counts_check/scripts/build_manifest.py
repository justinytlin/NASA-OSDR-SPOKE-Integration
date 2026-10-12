#!/usr/bin/env python3
"""Build the run manifest for the OSDR-wide tissue check from the survey cache.

One manifest entry = one (study, tissue) set of control samples:
  * spaceflight studies: samples whose Spaceflight factor is exactly "Ground Control"
  * ground-analog studies (hindlimb unloading, radiation, diet): samples at the
    baseline level of every condition factor (normally loaded, sham / non-irradiated,
    not injected, non-treated, normal diet)  -> control_type "analog control"
Genotype is restricted to wild type / control. Factors that are not conditions
(age, sex, strain, duration, dissection, mission, ...) are pooled.
Skipped with a reason: cell cultures, tumours, studies with no control group
(sample-preservation / methods studies), fewer than 2 control samples.

Writes data/_catalog/manifest.json and manifest.tsv (+ manifest_skipped.tsv).
"""
import json, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "_catalog"
ALREADY = {"OSD-564", "OSD-244", "OSD-245", "OSD-246", "OSD-288", "OSD-289",       # first set (kept as run)
           "OSD-137", "OSD-161", "OSD-162", "OSD-163", "OSD-168", "OSD-194", "OSD-270", "OSD-352"}   # RR-3 set

BASE = re.compile(r"^ground control$|normal(ly)? load(ed|ing) control|hindlimb loaded control|non-irradiated|"
                  r"sham[- ]irradiat|^not .*injected$|^wild type$|^control$|^non-treated$|^normal diet$|"
                  r"control diet containing 1% sucrose", re.I)
TREAT = re.compile(r"space flight|unloaded|radiation|proton|neutron|injected|ko$|null|mko|corticosterone|western diet|cape|"
                   r"basal control|vivarium control|cohort control|rerun|minimized mechanical stress|clinostat|nanoceria", re.I)
# tissue string -> (label, Anatomy family in families.py)
TISSUES = [
    (r"cell|cultur", None, None), (r"tumou?r", None, None), (r"plasma", None, None),
    (r"hippocamp", "hippocampus", "hippocampus"), (r"cerebell", "cerebellum", "cerebellum"),
    (r"cerebr(um|al)", "cerebrum", "cerebrum"), (r"brain", "brain", "brain"),
    (r"optic nerve", "optic nerve", "optic nerve"), (r"retina", "retina", "retina"), (r"\beye", "eye", "eye"),
    (r"liver", "liver", "liver"), (r"adrenal", "adrenal gland", "adrenal gland"), (r"kidney", "kidney", "kidney"),
    (r"spleen", "spleen", "spleen"), (r"thymus", "thymus", "thymus"), (r"bone marrow", "bone marrow", "bone marrow"),
    (r"extensor digitorum", "extensor digitorum longus", "skeletal muscle"), (r"gastrocnemius", "gastrocnemius", "skeletal muscle"),
    (r"quadriceps", "quadriceps", "skeletal muscle"), (r"soleus", "soleus", "skeletal muscle"),
    (r"tibialis", "tibialis anterior", "skeletal muscle"),
    (r"skin", "skin", "skin"), (r"colon", "colon", "colon"), (r"c(a)?ecum", "cecum", "cecum"), (r"lung", "lung", "lung"),
    (r"ventricle|heart", "heart", "heart"), (r"mammary", "mammary gland", "mammary gland"),
    (r"temporal bone", "temporal bone", "temporal bone"), (r"mandib", "mandible", "mandible"),
    (r"cortical bone|femur|tibia|\bbone", "bone", "bone tissue"),
    (r"brown adipose", "brown adipose tissue", "adipose tissue"), (r"white adipose", "white adipose tissue", "adipose tissue"),
    (r"adipose", "adipose tissue", "adipose tissue"), (r"blood", "blood", "whole blood"), (r"placenta", "placenta", "placenta"),
]


def norm_tissue(t):
    for rx, label, fam in TISSUES:
        if re.search(rx, str(t), re.I):
            return label, fam
    return "?", "?"


def main():
    survey = pd.read_csv(CAT / "survey.tsv", sep="\t").fillna("")
    entries, skipped = [], []
    for _, st in survey.iterrows():
        acc = st.accession
        if not st.counts_files:
            skipped.append((acc, st.material, "no GeneLab RSEM counts file")); continue
        if acc in ALREADY:
            continue
        s_files = list((CAT / "isa" / acc).rglob("s_*.txt"))
        df = pd.read_csv(s_files[0], sep="\t", dtype=str)
        df = df.loc[:, ~df.columns.duplicated()]
        cols = list(pd.read_csv(ROOT / "data" / st.counts_files, nrows=0).columns[1:])
        df = df[df["Sample Name"].isin(cols)].copy()
        if df.empty:
            skipped.append((acc, st.material, "ISA sample names do not match counts columns")); continue
        fcols = [c for c in df.columns if c.startswith("Factor Value[")]
        tcol = next((c for c in fcols if re.fullmatch(r"Factor Value\[Tissue\]", c, re.I)), None) or "Characteristics[Material Type]"
        keep = pd.Series(True, index=df.index); cond = []
        for c in fcols:
            if c == tcol:
                continue
            lv = df[c].fillna("").astype(str)
            is_cond = lv.map(lambda x: bool(BASE.search(x) or TREAT.search(x))).any()
            if is_cond:
                cond.append(c[len("Factor Value["):-1])
                keep &= lv.map(lambda x: bool(BASE.search(x)))
        if not cond:
            skipped.append((acc, st.material, "no condition factor (sample-preservation / methods study)")); continue
        ctrl = df[keep]
        ctype = "ground control" if any(re.fullmatch(r"Factor Value\[Spaceflight\]", c, re.I) for c in fcols) else "analog control"
        lab = ctrl[tcol].map(norm_tissue)
        ctrl = ctrl.assign(_label=[x[0] for x in lab], _family=[x[1] for x in lab])
        groups = ctrl.groupby(["_label", "_family"], dropna=False)
        multi = ctrl._label.nunique() > 1
        if ctrl.empty:
            skipped.append((acc, st.material, "no samples at the baseline level of " + ", ".join(cond))); continue
        for (label, fam), g in groups:
            if label is None or fam is None or pd.isna(label):
                skipped.append((acc, "; ".join(sorted(set(g[tcol]))), "cell culture / tumour / plasma")); continue
            if label == "?":
                skipped.append((acc, "; ".join(sorted(set(g[tcol].astype(str)))), "tissue not mapped")); continue
            if len(g) < 2:
                skipped.append((acc, label, f"only {len(g)} control sample")); continue
            name = f"{acc}_{label.replace(' ', '-')}" if multi else acc
            entries.append(dict(name=name, accession=acc, payload=st.payload, tissue=label, family=fam, control_type=ctype,
                                condition_factors=", ".join(cond), file=st.counts_files, n_controls=len(g),
                                raw_tissue="; ".join(sorted(set(g[tcol].astype(str)))), samples=list(g["Sample Name"])))
    (CAT / "manifest.json").write_text(json.dumps(entries, indent=1))
    m = pd.DataFrame(entries)
    m.drop(columns="samples").to_csv(CAT / "manifest.tsv", sep="\t", index=False)
    pd.DataFrame(skipped, columns=["accession", "tissue", "reason"]).to_csv(CAT / "manifest_skipped.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(m[["name", "payload", "tissue", "family", "control_type", "n_controls", "raw_tissue"]].to_string(index=False, max_colwidth=34))
    print(f"\n{len(m)} entries from {m.accession.nunique()} studies; {m.n_controls.sum()} control samples; "
          f"{(m.control_type == 'ground control').sum()} ground control / {(m.control_type == 'analog control').sum()} analog control")
    print("\ntissues:", m.groupby('tissue').size().sort_values(ascending=False).to_dict())
    print("\nskipped:"); print(pd.DataFrame(skipped, columns=["accession", "tissue", "reason"]).to_string(index=False, max_colwidth=70))


if __name__ == "__main__":
    main()

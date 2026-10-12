#!/usr/bin/env python3
"""Survey NASA OSDR for mouse bulk RNA-seq studies usable by the GC tissue check.

For every Mus musculus study with an RNA-seq assay: cache the file listing and
the ISA metadata, find the GeneLab RSEM unnormalized counts file, and tabulate
the study factors, their levels, and the organism part of every sample.
Writes data/_catalog/{files_<OSD>.json, isa/<OSD>/...} and
data/_catalog/survey.tsv (one row per study) + survey_samples.tsv (one row per sample).

Usage: survey_osdr.py [--refresh]
"""
import argparse, io, json, re, sys, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import urllib.parse, urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CAT = ROOT / "data" / "_catalog"
API = "https://osdr.nasa.gov"


def get(url, binary=False, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                b = r.read()
            return b if binary else b.decode("utf8", "replace")
        except Exception as e:  # noqa: BLE001
            err = e
    raise err


def study_list():
    q = urllib.parse.urlencode({"term": "", "type": "cgene", "size": 1000, "ffield": "organism", "fvalue": "Mus musculus"})
    d = json.loads(get(f"{API}/osdr/data/search?{q}"))
    out = []
    for h in d["hits"]["hits"]:
        s = h["_source"]
        if "RNA Sequencing (RNA-Seq)" in str(s.get("Study Assay Technology Type")):
            out.append(dict(accession=s["Accession"], payload=str(s.get("Project Identifier") or ""),
                            title=str(s.get("Study Title") or ""), material=re.sub(r"\s{2,}", "; ", str(s.get("Material Type") or ""))))
    return sorted(out, key=lambda r: int(r["accession"].split("-")[1]))


def fetch(st, refresh):
    acc = st["accession"]; n = acc.split("-")[1]
    fj = CAT / f"files_{acc}.json"
    if refresh or not fj.exists():
        fj.write_text(get(f"{API}/osdr/data/osd/files/{n}"))
    files = json.loads(fj.read_text())["studies"][acc]["study_files"]
    names = [f["file_name"] for f in files]
    counts = [x for x in names if re.search(r"rna_seq_RSEM_Unnormalized_Counts(_GLbulkRNAseq)?\.csv$", x) and "rRNArm" not in x]
    isa = [x for x in names if x.endswith("-ISA.zip")]
    idir = CAT / "isa" / acc
    if isa and (refresh or not idir.exists()):
        b = get(f"{API}/geode-py/ws/studies/{acc}/download?source=datamanager&file={urllib.parse.quote(isa[0])}", binary=True)
        idir.mkdir(parents=True, exist_ok=True)
        zipfile.ZipFile(io.BytesIO(b)).extractall(idir)
    st.update(counts_files="; ".join(counts), n_files=len(names))
    return st


def sample_table(acc):
    idir = CAT / "isa" / acc
    s_files = list(idir.rglob("s_*.txt"))
    if not s_files:
        return None
    df = pd.read_csv(s_files[0], sep="\t", dtype=str)
    keep = {"Sample Name": "sample"}
    for c in df.columns:
        lc = c.lower()
        if c.startswith("Factor Value["):
            keep[c] = "F:" + c[len("Factor Value["):].rstrip("]")
        elif "organism part" in lc or lc.startswith("characteristics[material type"):
            keep[c] = "tissue" if "organism part" in lc else "material_type"
        elif lc.startswith("characteristics[strain") or lc.startswith("characteristics[genotype") or lc.startswith("characteristics[sex"):
            keep[c] = c[len("Characteristics["):].rstrip("]").lower()
    out = df[list(keep)].rename(columns=keep)
    out = out.loc[:, ~out.columns.duplicated()]
    out.insert(0, "accession", acc)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--refresh", action="store_true"); a = ap.parse_args()
    CAT.mkdir(parents=True, exist_ok=True)
    studies = study_list()
    print(f"{len(studies)} mouse studies with an RNA-seq assay", flush=True)
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda s: _safe(fetch, s, a.refresh), studies))
    rows, samples = [], []
    for st in res:
        acc = st["accession"]
        t = sample_table(acc) if "error" not in st else None
        if t is not None:
            samples.append(t)
            fcols = [c for c in t.columns if c.startswith("F:")]
            st["n_samples"] = len(t)
            st["tissues"] = "; ".join(sorted(set(t["tissue"].dropna()))) if "tissue" in t else ""
            st["factors"] = " || ".join(f"{c[2:]}: " + ", ".join(f"{k} ({v})" for k, v in t[c].value_counts(dropna=False).items()) for c in fcols)
        rows.append(st)
    pd.DataFrame(rows).to_csv(CAT / "survey.tsv", sep="\t", index=False)
    pd.concat(samples, ignore_index=True).to_csv(CAT / "survey_samples.tsv", sep="\t", index=False)
    df = pd.DataFrame(rows)
    print("with GeneLab RSEM unnormalized counts:", int((df.counts_files.fillna("") != "").sum()))
    print("wrote", CAT / "survey.tsv")


def _safe(fn, st, refresh):
    try:
        return fn(st, refresh)
    except Exception as e:  # noqa: BLE001
        st["error"] = repr(e)[:200]; st.setdefault("counts_files", ""); return st


if __name__ == "__main__":
    main()

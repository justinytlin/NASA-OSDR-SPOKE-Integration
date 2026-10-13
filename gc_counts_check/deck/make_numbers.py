"""Collect the edge-set comparison numbers for build.js -> numbers.json (run after pilot_stage2.py for all four tags)."""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results"
tags = ["2025", "2025noAeG", "2025nohk", "2025upAuG"]
T = {t: pd.read_csv(R / f"pilot_stage2_{t}.tsv", sep="\t").set_index("tissue") for t in tags}
tissues = list(T["2025"].index)
def cell(t, tissue):
    r = T[t].loc[tissue, f"rank_{t}"]; v = str(T[t].loc[tissue, f"verdict_{t}"])
    if pd.isna(r): return "no node"
    return f"{int(r)}" + (" (proxy)" if "proxy" in v else "")
rows = [[tissue] + [cell(t, tissue) for t in tags] for tissue in tissues]
nl = {t: int(T[t][f"n_linked_{t}"].dropna().iloc[0]) for t in tags}
# connectivity correlation per build: link weight (mean PSEV at the node, from _anatomy_psev_<tag>.npz) vs own-node rank
def spear(t):
    A = np.load(ROOT / "psev_out" / f"_anatomy_psev_{t}.npz"); V = A["values"]; link = (V > 1e-4).mean(axis=0) > 0.5
    build = {"2025": "gene_psev_spoke2025-03", "2025noAeG": "gene_psev_spoke2025-03_noAeG", "2025nohk": "gene_psev_spoke2025-03", "2025upAuG": "gene_psev_spoke2025-03_upAuG_only"}[t]
    ni = pd.read_csv(f"/Volumes/justinytlin/SPOKE_psev/{build}/node_info.tsv", sep="\t", dtype={"Node": str}); an = ni[ni.Node_Type == "Anatomy"].reset_index(drop=True)
    lw = pd.Series((-V.mean(axis=0)[link]).argsort().argsort() + 1, index=an.Node_Name.astype(str).values[link])
    xs, ys = [], []
    for tissue in tissues:
        tg = T[t].loc[tissue, f"targets_{t}"]; r = T[t].loc[tissue, f"rank_{t}"]; v = str(T[t].loc[tissue, f"verdict_{t}"])
        if pd.isna(r) or "proxy" in v or not isinstance(tg, str): continue
        name = re.sub(r" \(rank \d+\)", "", tg.split(";")[0]).strip()
        if name in lw.index: xs.append(lw[name]); ys.append(r)
    return f"{spearmanr(xs, ys).correlation:+.2f}".replace("+", "") if len(xs) > 3 else "n/a"
sp = {t: spear(t) for t in tags}
lf = int(T["2025"].loc["liver", "rank_2025"]); ordn = lambda n: f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"
out = dict(stage2_rows=rows, n_linked_full=nl["2025"], n_linked_noAeG=nl["2025noAeG"], n_linked_nohk=nl["2025nohk"], n_linked_up=nl["2025upAuG"],
           spearman_full=sp["2025"].replace("-", "−"), spearman_noAeG=sp["2025noAeG"].replace("-", "−"), spearman_nohk=sp["2025nohk"].replace("-", "−"), spearman_up=sp["2025upAuG"].replace("-", "−"),
           liver_full=lf, liver_full_ordinal=ordn(lf), spleen_full=cell("2025", "spleen"), blood_full=cell("2025", "blood"), kidney_full=cell("2025", "kidney"))
json.dump(out, open(Path(__file__).resolve().parent / "numbers.json", "w"), indent=1); print(json.dumps(out, indent=1))

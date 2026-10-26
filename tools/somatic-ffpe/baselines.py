"""Author baselines emulating common (incomplete) pipelines on cached fragment counts.

python3 baselines.py PUBLIC PRIVATE  -> precision/recall and naive purity/ploidy/TiN per baseline
"""
import csv, json, sys
from pathlib import Path
import numpy as np
import pysam
sys.path.insert(0, "solution")
from call_somatic import count_library, B2I
pub, priv = Path(sys.argv[1]), Path(sys.argv[2])
fa = pysam.FastaFile(str(pub / "reference.fa"))
names = list(fa.references); lens = [fa.get_reference_length(c) for c in names]
starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(int); offsets = dict(zip(names, starts)); L = int(sum(lens))
cache = priv / "counts.npz"
if not cache.exists():
    t, td = count_library(str(pub / "tumor.cram"), str(pub / "reference.fa"), offsets, L)
    n, nd = count_library(str(pub / "normal.cram"), str(pub / "reference.fa"), offsets, L)
    np.savez_compressed(cache, t=t, n=n, td=td, nd=nd)
d = np.load(cache); t, n, td, nd = d["t"], d["n"], d["td"], d["nd"]
ref = np.concatenate([B2I[np.frombuffer(fa.fetch(c).encode(), dtype=np.uint8)] for c in names])
P = np.arange(L)
truth = set()
for r in csv.DictReader(open(priv / "truth_somatic.tsv"), delimiter="\t"):
    truth.add((offsets[r["chrom"]] + int(r["pos"]) - 1, r["alt"]))
tp_ = json.load(open(priv / "truth_params.json"))
ta = t.sum(0); na = n.sum(0); tdp = ta.sum(1); ndp = na.sum(1)
x = ta.copy(); x[P, np.minimum(ref, 3)] = -1; alt = x.argmax(1)
at = ta[P, alt]; an = na[P, alt]; vt = at / np.maximum(tdp, 1); vn = an / np.maximum(ndp, 1)
prone = ((ref == 1) & (alt == 3)) | ((ref == 2) & (alt == 0))
bad_o = np.where((ref == 1) & (alt == 3), 0, 1)
clean_alt = np.where(prone, t[1 - bad_o, P, alt], at)
base = (ref < 4) & (at >= 4) & (vt >= 0.03)
orient_ok = ~prone | (clean_alt >= 2)
def report(name, mask, tin=None):
    calls = {(int(p), "ACGT"[alt[p]]) for p in np.flatnonzero(mask)}
    tp = len(calls & truth); prec = tp / max(len(calls), 1); rec = tp / len(truth)
    hi = np.flatnonzero(mask & (vt > 0.15))
    pur = min(1, 2 * float(np.median(vt[hi]))) if len(hi) else float("nan")
    print(f"{name:34s} calls={len(calls):5d} P={prec:.3f} R={rec:.3f} naive_purity(2xmedian)={pur:.3f}" + (f" tin={tin:.3f}" if tin is not None else ""))
report("B1 strict normal, no FFPE", base & (an <= 1))
report("B2 orientation + strict N", base & orient_ok & (an <= 1))
report("B3 orientation + germline-only N<0.2", base & orient_ok & (vn < 0.2))
report("B3b orientation + germline-only N<0.3", base & orient_ok & (vn < 0.3))
cl = base & orient_ok & (vn < 0.3) & (vt > 0.15)
rho = float(np.median(vn[cl] / vt[cl]))
report("B4 TiN ratio, no FFPE", base & (vn < 0.3) & (vn <= 0.6 * vt), rho * 0.62)
report("B5 TiN ratio, drop lowVAF C>T", base & (vn < 0.3) & (vn <= 0.6 * vt) & ~(prone & (vt < 0.15)), rho * 0.62)
report("B6 TiN ratio + orientation", base & orient_ok & (vn < 0.3) & (vn <= 0.6 * vt), rho * 0.62)
# B8: emulation of round-1 trial-1 filters: clean-orientation evidence pooled over tumour and
# normal, germline normal VAF < 0.35, removal of mate-discordant sites, no N/T ratio model
clean_n = na[P, alt]
tot = clean_alt + clean_n
disc = td[P, alt] + nd[P, alt]
cand8 = (ref < 4) & (vn < 0.35) & (disc <= 1) & np.where(prone, tot >= 3, tot >= 4) & (at >= 2)
report("B8 trial-1 emulation (no ratio)", cand8)
print("truth purity", tp_["tumor_purity"], "ploidy", tp_["tumor_ploidy"], "tin", tp_["normal_tumor_fraction"], "| median N/T ratio", round(rho, 3))

import csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "solution")
from call_somatic import count_library
import pysam
pub, priv, res = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cache = priv / "counts.npz"
if not cache.exists():
    L = pysam.FastaFile(str(pub / "reference.fa")).get_reference_length("chrS")
    tt, td = count_library(str(pub / "tumor.cram"), str(pub / "reference.fa"), L)
    nn, nd = count_library(str(pub / "normal.cram"), str(pub / "reference.fa"), L)
    np.savez_compressed(cache, t=tt, n=nn, td=td, nd=nd)
d = np.load(cache); t, n = d["t"], d["n"]
bi = {b: i for i, b in enumerate("ACGT")}
truth = {int(r["pos"]) for r in csv.DictReader(open(priv / "truth_somatic.tsv"), delimiter="\t")}
rec = {int(r["pos"]): r for r in csv.DictReader(open(priv / "truth_recurrent.tsv"), delimiter="\t")}
for r in csv.DictReader(open(res / "somatic_snvs.tsv"), delimiter="\t"):
    p = int(r["pos"]); a = bi[r["alt"]]
    if p in truth: continue
    tag = "REC rate=" + rec[p]["rate"] if p in rec else "other"
    print(p, r["ref"] + ">" + r["alt"], "T", t[0, p-1, a], t[0, p-1].sum(), t[1, p-1, a], t[1, p-1].sum(), "N", n[:, p-1, a].sum(), n[:, p-1].sum(), tag)

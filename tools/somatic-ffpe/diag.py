"""Author diagnostics: cache fragment counts and print features at truth/artifact sites."""
import csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "solution")
from call_somatic import count_library, B2I
import pysam
pub, priv = Path(sys.argv[1]), Path(sys.argv[2])
cache = priv / "counts.npz"
fa = pysam.FastaFile(str(pub / "reference.fa")); L = fa.get_reference_length("chrS")
if not cache.exists():
    t = count_library(str(pub / "tumor.cram"), str(pub / "reference.fa"), L)
    n = count_library(str(pub / "normal.cram"), str(pub / "reference.fa"), L)
    np.savez_compressed(cache, t=t, n=n)
d = np.load(cache); t, n = d["t"], d["n"]
bi = {b: i for i, b in enumerate("ACGT")}
print("clone pos ref>alt | tF1R2 alt/dep tF2R1 alt/dep | normal alt/dep")
for r in csv.DictReader(open(priv / "truth_somatic.tsv"), delimiter="\t"):
    p = int(r["pos"]) - 1; a = bi[r["alt"]]
    row = (r["clone"], p + 1, r["ref"] + ">" + r["alt"], t[0, p, a], t[0, p].sum(), t[1, p, a], t[1, p].sum(), n[:, p, a].sum(), n[:, p].sum())
    print(*row)

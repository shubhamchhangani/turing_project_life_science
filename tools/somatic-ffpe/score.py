"""Author scoring helper: python3 score.py PRIVATE_DIR RESULTS_DIR"""
import csv, json, sys
from pathlib import Path
priv, res = Path(sys.argv[1]), Path(sys.argv[2])
truth = {}
for r in csv.DictReader(open(priv / "truth_somatic.tsv"), delimiter="\t"):
    truth[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = r["clone"]
rec = {(r["chrom"], int(r["pos"])) for r in csv.DictReader(open(priv / "truth_recurrent.tsv"), delimiter="\t")}
germ = {(r["chrom"], int(r["pos"])) for r in csv.DictReader(open(priv / "truth_germline.tsv"), delimiter="\t")}
calls = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(open(res / "somatic_snvs.tsv"), delimiter="\t")}
tp = calls & truth.keys()
fp = calls - truth.keys()
P = len(tp) / max(len(calls), 1); R = len(tp) / len(truth)
by = {c: (sum(1 for k in tp if truth[k] == c), sum(1 for v in truth.values() if v == c)) for c in ("clonal", "A", "B")}
fpk = {"recurrent": sum((k[0], k[1]) in rec for k in fp), "germline": sum((k[0], k[1]) in germ for k in fp)}
fpk["ffpe_like"] = sum(1 for k in fp if (k[2], k[3]) in (("C", "T"), ("G", "A")) and (k[0], k[1]) not in rec and (k[0], k[1]) not in germ)
fpk["other"] = len(fp) - sum(fpk.values())
tp_ = json.load(open(priv / "truth_params.json"))
est = json.load(open(res / "sample_estimates.json"))
print(f"calls={len(calls)} TP={len(tp)} P={P:.3f} R={R:.3f} by_clone={by} FP={fpk}")
print(f"ploidy {est.get('tumor_ploidy',float('nan')):.3f} vs {tp_.get('tumor_ploidy',2.0)}; purity {est['tumor_purity']:.4f} vs {tp_['tumor_purity']} (rel {est['tumor_purity']/tp_['tumor_purity']-1:+.3f}); "
      f"tin {est['normal_tumor_fraction']:.4f} vs {tp_['normal_tumor_fraction']} (rel {est['normal_tumor_fraction']/tp_['normal_tumor_fraction']-1:+.3f})")

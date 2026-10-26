"""Authoring tool: find the first seed whose estimates are well separated from the
0.05 decision threshold and from zero (so relative tolerance is meaningful), and
which hits the intended pattern A: CT only, B: SC only, C: both, D: neither.

Usage: python pick_seed.py <start_seed> <work_dir>
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
PATTERN = {"Locus_A": ("yes", "no"), "Locus_B": ("no", "yes"), "Locus_C": ("yes", "yes"), "Locus_D": ("no", "no")}
MIN_GAP = 0.02      # |F - 0.05| for F_CT and F_SC
MIN_ABS = 0.005     # |F| for every statistic


def main():
    seed = int(sys.argv[1])
    work = Path(sys.argv[2]); work.mkdir(parents=True, exist_ok=True)
    csv_path = work / "genotype_counts.csv"
    while True:
        subprocess.run([sys.executable, HERE / "make_counts.py", str(seed), csv_path], check=True)
        res = json.loads(subprocess.run([sys.executable, HERE / "bruteforce_fstats.py", csv_path],
                                        check=True, capture_output=True, text=True).stdout)
        ok = True
        for loc, f in res.items():
            pat = ("yes" if f["F_CT"] > 0.05 else "no", "yes" if f["F_SC"] > 0.05 else "no")
            ok &= pat == PATTERN[loc]
            ok &= all(abs(f[k] - 0.05) >= MIN_GAP for k in ("F_CT", "F_SC"))
            ok &= all(abs(v) >= MIN_ABS for v in f.values())
        if ok:
            print(f"seed {seed}")
            print(json.dumps(res, indent=2))
            return
        seed += 1


if __name__ == "__main__":
    main()

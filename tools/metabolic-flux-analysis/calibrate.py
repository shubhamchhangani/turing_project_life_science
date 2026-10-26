"""Authoring tool: Monte-Carlo calibration of grading tolerances.

Generates N independent noisy datasets from the true fluxes (seeds 1000..1000+N-1),
fits each with the reference solution, and reports the spread of the estimated net
fluxes around the truth. Grading tolerances are set from this spread, so a correct
analysis of any noise realization passes.

Usage: python calibrate.py <n_datasets> <n_workers> <out_json>
"""
import csv
import json
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
FIT = ROOT / "tasks/life-sciences/biology/metabolic-flux-analysis/solution/mfa_fit.py"
PY = sys.executable
GRADED = ["G6PDH", "PDH", "PC", "ME", "CS", "IDH", "GLUDH", "LDH", "PGI", "MDH"]


def one(seed):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        subprocess.run([PY, HERE / "make_dataset.py", str(seed), tmp / "data", tmp / "truth.json"],
                       check=True, capture_output=True)
        subprocess.run([PY, FIT, tmp / "data", tmp / "out.csv", "4", str(seed)],
                       check=True, capture_output=True)
        with open(tmp / "out.csv") as fh:
            rows = {r["reaction_id"]: r for r in csv.DictReader(fh)}
    return seed, {k: float(rows[k]["net_flux"]) for k in rows}, \
        {k: float(rows[k]["exchange_flux"]) for k in rows if rows[k]["exchange_flux"]}


def main():
    n, workers, out = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
    sys.path.insert(0, str(HERE))
    import true_fluxes as T
    with ProcessPoolExecutor(workers) as ex:
        results = list(ex.map(one, range(1000, 1000 + n)))
    summary = {}
    for r in GRADED:
        est = np.array([res[1][r] for res in results])
        err = est - T.NET[r]
        summary[r] = {"true": T.NET[r], "mean": est.mean(), "sd": est.std(ddof=1),
                      "max_abs_err": np.abs(err).max(), "max_rel_err": (np.abs(err) / T.NET[r]).max()}
        print(f"{r:6s} true {T.NET[r]:7.2f}  mean {est.mean():7.2f}  sd {est.std(ddof=1):6.2f}  "
              f"max|err| {np.abs(err).max():6.2f}  max rel {summary[r]['max_rel_err']:.3f}")
    for r in sorted(results[0][2]):
        est = np.array([res[2][r] for res in results])
        print(f"exch {r:6s} true {T.EXCH.get(r, 0):7.2f}  mean {est.mean():7.2f}  sd {est.std(ddof=1):6.2f}")
    out.write_text(json.dumps({"summary": summary, "runs": results}, indent=1))


if __name__ == "__main__":
    main()

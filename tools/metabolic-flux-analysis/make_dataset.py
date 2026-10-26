"""Authoring tool: write the task inputs (agent-visible) and the truth file (verifier only).

Usage: python make_dataset.py <seed> <env_data_dir> <truth_json>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

import mfa_model as M
import true_fluxes as T

MID_SD = 0.003
EXTRACELLULAR = {   # reaction id: relative SD of the measurement (glucose uptake is the exact reference)
    "LACout": 0.05, "GLNup": 0.05, "GLUout": 0.10, "ALAout": 0.10, "SERout": 0.10,
    "ASPout": 0.10, "P5Pout": 0.10, "ACCOAout": 0.10,
}


def equation(rid):
    subs, prods, rev = M.REACTIONS[rid]
    lhs = " + ".join(m for m, _ in subs)
    if rid == "GLCup":
        lhs = "Glc.ext"
    if rid == "GLNup":
        lhs = "Gln.ext"
    rhs = " + ".join(m for m, _ in prods) if prods else subs[0][0] + ".ext"
    arrow = " <=> " if rev else " -> "
    return lhs + arrow + rhs


def formula_str(form):
    return "".join(el + (str(form[el]) if form[el] > 1 else "") for el in ("C", "H", "N", "O", "Si") if form.get(el))


def main():
    seed, data_dir, truth_path = int(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    data_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    with open(data_dir / "network.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "equation", "reversible"])
        for rid, (_s, _p, rev) in M.REACTIONS.items():
            w.writerow([rid, equation(rid), "yes" if rev else "no"])

    with open(data_dir / "fragments.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fragment", "metabolite", "backbone_carbons", "ion_formula"])
        for f, (met, carbons, form) in M.FRAGMENTS.items():
            w.writerow([f, met, ";".join(map(str, carbons)), formula_str(form)])

    with open(data_dir / "isotope_abundances.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "mass_shift", "abundance"])
        for el, iso in M.ISOTOPES.items():
            for shift, ab in iso:
                w.writerow([el, shift, ab])

    with open(data_dir / "mid_measurements.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["experiment", "fragment", "isotopologue", "fraction", "sd"])
        for exp in M.EXPERIMENTS:
            dists = M.simulate(T.NET, T.EXCH, M.sources_for(exp))
            for frag in M.FRAGMENTS:
                mid = M.raw_mid(dists, frag)
                noisy = mid + rng.normal(0.0, MID_SD, mid.size)
                for k, val in enumerate(noisy):
                    w.writerow([exp, frag, f"M{k}", f"{val:.4f}", f"{MID_SD:.4f}"])

    with open(data_dir / "extracellular_rates.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction_id", "rate", "sd"])
        w.writerow(["GLCup", "100.0", "0"])
        for rid, rel in EXTRACELLULAR.items():
            sd = rel * T.NET[rid]
            w.writerow([rid, f"{T.NET[rid] + rng.normal(0.0, sd):.2f}", f"{sd:.2f}"])

    truth = {"net": T.NET, "exchange": T.EXCH, "seed": seed}
    truth_path.write_text(json.dumps(truth, indent=2) + "\n")
    print("wrote", data_dir, truth_path)


if __name__ == "__main__":
    main()

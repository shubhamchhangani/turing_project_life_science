"""Authoring tool: negative controls. Fits the task data with plausible-but-wrong analyses.

Each variant changes exactly one thing in the reference solution (by monkeypatching
or by transforming the input data the way a careless analyst would) and writes a
fluxes.csv that is then graded by the verifier.

Usage: python wrong_methods.py <data_dir> <out_dir> [variant ...]
"""
import csv
import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIT = ROOT / "tasks/life-sciences/biology/metabolic-flux-analysis/solution/mfa_fit.py"


def load_fit(data_dir, out_csv):
    sys.argv = ["mfa_fit.py", str(data_dir), str(out_csv), "4", "11"]
    spec = importlib.util.spec_from_file_location("mfa_fit_variant", FIT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def copy_data(src, dst):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def rewrite_csv(path, fn, delimiter=","):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter=delimiter))
        fields = list(rows[0])
    rows = fn(rows)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter=delimiter)
        w.writeheader()
        w.writerows(rows)


def no_derivatization_correction(d):
    """Treat raw MIDs as backbone MIDs: keep M+0..M+n renormalized, ignore other atoms."""
    nb = {}
    def frag_fn(rows):
        for r in rows:
            nb[r["fragment"]] = len(r["backbone_carbons"].split(";"))
            r["ion_formula"] = f"C{nb[r['fragment']]}"
        return rows
    rewrite_csv(d / "fragments.csv", frag_fn)

    def mid_fn(rows):
        keep = [r for r in rows if int(r["isotopologue"][1:]) <= nb[r["fragment"]]]
        tot = {}
        for r in keep:
            tot[(r["experiment"], r["fragment"])] = tot.get((r["experiment"], r["fragment"]), 0) + float(r["fraction"])
        for r in keep:
            r["fraction"] = f"{float(r['fraction']) / tot[(r['experiment'], r['fragment'])]:.4f}"
        return keep
    rewrite_csv(d / "mid_measurements.csv", mid_fn)


def single_experiment(name):
    def fn(d):
        rewrite_csv(d / "mid_measurements.csv", lambda rows: [r for r in rows if r["experiment"] == name])
    return fn


def idh_irreversible(d):
    def fn(rows):
        for r in rows:
            if r["reaction_id"] == "IDH":
                r["reversible"] = "no"
        return rows
    rewrite_csv(d / "network.tsv", fn, delimiter="\t")


def pure_12c(mod):
    def sources():
        z = lambda n: mod.positional([0.0] * n)
        return {"glc12": {"Glc": mod.positional([0.99, 0.99, 0, 0, 0, 0]), "Gln": z(5), "CO2": z(1)},
                "gln5": {"Glc": z(6), "Gln": mod.positional([0.99] * 5), "CO2": z(1)}}
    mod.tracer_sources = sources


def perfect_purity(mod):
    def sources():
        nat = mod.natural
        return {"glc12": {"Glc": mod.positional([1.0, 1.0] + [0.0107] * 4), "Gln": nat(5), "CO2": nat(1)},
                "gln5": {"Glc": nat(6), "Gln": mod.positional([1.0] * 5), "CO2": nat(1)}}
    mod.tracer_sources = sources


def no_symmetry(mod):
    mod.SYMMETRIC = set()


def naive_cs(mod):
    mod.ATOMS["CS"] = ([("AcCoA", "ab"), ("OAA", "cdef")], [("Cit", "abcdef")])


VARIANTS = {
    "no_derivatization_correction": ("data", no_derivatization_correction),
    "unlabeled_carbon_as_pure_12C": ("mod", pure_12c),
    "no_succinate_fumarate_symmetry": ("mod", no_symmetry),
    "glucose_tracer_only": ("data", single_experiment("glc12")),
    "glutamine_tracer_only": ("data", single_experiment("gln5")),
    "naive_citrate_synthase_map": ("mod", naive_cs),
    "idh_irreversible": ("data", idh_irreversible),
    "tracer_purity_ignored": ("mod", perfect_purity),
}


def run(variant, data_dir, out_dir):
    kind, fn = VARIANTS[variant]
    vdir = out_dir / variant
    vdir.mkdir(parents=True, exist_ok=True)
    d = vdir / "data"
    copy_data(data_dir, d)
    if kind == "data":
        fn(d)
    mod = load_fit(d, vdir / "results" / "fluxes.csv")
    if kind == "mod":
        fn(mod)
    mod.main()


if __name__ == "__main__":
    data_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    names = sys.argv[3:] or list(VARIANTS)
    for v in names:
        print("=====", v, flush=True)
        run(v, data_dir, out_dir)

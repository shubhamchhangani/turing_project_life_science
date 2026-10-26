"""Authoring tool: cross-validate mfa_model.simulate against FreeFlux (Wu et al. 2023).

Runs under the FreeFlux venv (Python 3.10). Builds the same network in FreeFlux
format (EMU-based solver, completely independent code), simulates the backbone
MDVs of every measured fragment at the true fluxes for both tracer experiments
and prints them as JSON for comparison.

Usage: python crosscheck_freeflux.py <out_json>
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import importlib, os  # noqa: E402
T = importlib.import_module(os.environ.get("FLUXES", "true_fluxes"))

from freeflux import Model  # noqa: E402

# Pentose and triose pools are split into separately named pools joined by near-equilibrium
# isomerase reactions, because FreeFlux indexes reactants by name and cannot hold the same
# metabolite twice on one side of a reaction (P5P + P5P, F6P -> T3P + T3P).
# Carboxylations take CO2 from an unlabeled medium pool (CO2in); released CO2 leaves (CO2out).
MODEL_TSV = """#reaction_ID	substrate_IDs(atom)	product_IDs(atom)	reversibility
GLCup	Glc(abcdef)	G6P(abcdef)	0
PGI	G6P(abcdef)	F6P(abcdef)	1
PFK	F6P(abcdef)	DHAP(abc)+T3P(def)	0
TPI	DHAP(abc)	T3P(cba)	1
G6PDH	G6P(abcdef)	P5P(bcdef)+CO2out(a)	0
RPE	P5P(abcde)	X5P(abcde)	1
RPI	P5P(abcde)	R5P(abcde)	1
TK1	X5P(abcde)+R5P(fghij)	T3P(cde)+S7P(abfghij)	1
TA	S7P(abcdefg)+T3P(hij)	E4P(defg)+F6P(abchij)	1
TK2	X5P(abcde)+E4P(fghi)	T3P(cde)+F6P(abfghi)	1
GAPDH	T3P(abc)	PG3(abc)	0
PK	PG3(abc)	Pyr(abc)	0
LDH	Pyr(abc)	Lac(abc)	0
ALT	Pyr(abc)	Ala(abc)	0
PHGDH	PG3(abc)	Ser(abc)	0
PDH	Pyr(abc)	AcCoA(bc)+CO2out(a)	0
PC	Pyr(abc)+CO2in(d)	OAA(abcd)	0
CS	AcCoA(ab)+OAA(cdef)	Cit(fedbac)	0
IDHf	Cit(abcdef)	AKG(abcde)+CO2out(f)	0
IDHr	AKG(abcde)+CO2in(f)	Cit(abcdef)	0
AKGDH	AKG(abcde)	Suc(bcde,edcb)+CO2out(a)	0
SDH	Suc(abcd,dcba)	Fum(abcd,dcba)	0
FH	Fum(abcd,dcba)	Mal(abcd)	1
MDH	Mal(abcd)	OAA(abcd)	1
ME	Mal(abcd)	Pyr(abc)+CO2out(d)	0
GOT	OAA(abcd)	Asp(abcd)	0
GLNup	Gln(abcde)	Glu(abcde)	0
GLUDH	Glu(abcde)	AKG(abcde)	1
LACout	Lac(abc)	Lacex(abc)	0
ALAout	Ala(abc)	Alaex(abc)	0
SERout	Ser(abc)	Serex(abc)	0
ASPout	Asp(abcd)	Aspex(abcd)	0
GLUout	Glu(abcde)	Gluex(abcde)	0
P5Pout	P5P(abcde)	P5Pex(abcde)	0
ACCOAout	AcCoA(ab)	AcCoAex(ab)	0
"""

TARGETS = {"Lac": ["123"], "Ala": ["123", "23"], "Ser": ["123"], "Asp": ["1234", "234"],
           "Glu": ["12345", "2345"], "Fum": ["1234"], "Mal": ["1234"], "Cit": ["123456"]}


def directional():
    out = {}
    rev = {"PGI", "TK1", "TA", "TK2", "FH", "MDH", "GLUDH"}
    # isomerases/epimerase that FreeFlux needs instead of duplicate pools; near-equilibrium
    big = 1.0e5
    x5p_net = T.NET["TK1"] + T.NET["TK2"]
    r5p_net = T.NET["TK1"]
    out["RPE_f"], out["RPE_b"] = big + max(x5p_net, 0), big + max(-x5p_net, 0)
    out["RPI_f"], out["RPI_b"] = big + max(r5p_net, 0), big + max(-r5p_net, 0)
    out["TPI_f"], out["TPI_b"] = big + T.NET["PFK"], big
    for rid, v in T.NET.items():
        if rid == "IDH":
            x = T.EXCH["IDH"]
            out["IDHf"] = x + max(v, 0)
            out["IDHr"] = x + max(-v, 0)
        elif rid in rev:
            x = T.EXCH.get(rid, 0.0)
            out[rid + "_f"] = x + max(v, 0)
            out[rid + "_b"] = x + max(-v, 0)
        else:
            out[rid] = v
    return out


def run(experiment):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(MODEL_TSV)
        path = fh.name
    model = Model("crosscheck")
    model.read_from_file(path)
    with model.simulator("ss") as sim:
        sim.set_target_EMUs(TARGETS)
        if experiment == "glc12":
            sim.set_labeling_strategy("Glc", ["110000"], [1.0], [0.99])
        else:
            sim.set_labeling_strategy("Gln", ["11111"], [1.0], [0.99])
        for fid, v in directional().items():
            sim.set_flux(fid, v)
        sim.prepare()
        res = sim.simulate()
    out = {}
    for met, atoms_list in TARGETS.items():
        for atoms in atoms_list:
            emu = f"{met}_{atoms}"
            out[emu] = [float(x) for x in res._simulated_MDVs[emu].value]
    return out


if __name__ == "__main__":
    result = {e: run(e) for e in ("glc12", "gln5")}
    Path(sys.argv[1]).write_text(json.dumps(result, indent=1))
    print("ok")

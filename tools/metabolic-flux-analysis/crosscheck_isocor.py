"""Authoring tool: cross-validate the raw-MID measurement model against IsoCor 2.2.2
(Millard et al. 2019, Bioinformatics 35:4484).

For each fragment, the simulated raw MID (backbone MID convolved with the natural
isotope distribution of every non-backbone atom, truncated to M+0..M+(n+2) and
renormalized) is corrected with IsoCor's low-resolution corrector, treating the
non-backbone carbons as a derivative moiety and leaving the backbone carbons
uncorrected. The corrected vector must reproduce the simulated backbone MID.
"""
import numpy as np
from isocor.mscorrectors import LowResMetaboliteCorrector

import mfa_model as M
import true_fluxes as T

MASS = {"C": [12.0, 13.0033548378], "H": [1.0078250321, 2.0141017780],
        "N": [14.0030740052, 15.0001088984], "O": [15.9949146221, 16.99913150, 17.9991604],
        "Si": [27.9769265327, 28.97649472, 29.97377022]}
ISODATA = {el: {"abundance": [a for _, a in iso], "mass": MASS[el]} for el, iso in M.ISOTOPES.items()}


def formula_str(d):
    return "".join(f"{el}{n}" for el, n in d.items() if n > 0)


worst = 0.0
for exp in M.EXPERIMENTS:
    dists = M.simulate(T.NET, T.EXCH, M.sources_for(exp))
    for frag, (met, carbons, form) in M.FRAGMENTS.items():
        n = len(carbons)
        moiety = {"C": n, **{el: c for el, c in form.items() if el != "C"}}
        deriv = {"C": form["C"] - n}
        corr = LowResMetaboliteCorrector(formula_str(moiety), "13C", data_isotopes=ISODATA,
                                         derivative_formula=formula_str(deriv),
                                         correct_NA_tracer=False, tracer_purity=[0.0, 1.0])
        raw = M.raw_mid(dists, frag)
        _area, mid_iso, *_ = corr.correct(raw[: n + 1])   # IsoCor takes M+0..M+n
        mine = M.backbone_mid(dists[met], carbons)
        diff = np.max(np.abs(np.asarray(mid_iso) - mine))
        worst = max(worst, diff)
        print(f"{exp:6s} {frag:7s} max|IsoCor - model| = {diff:.1e}")
print("WORST", worst)

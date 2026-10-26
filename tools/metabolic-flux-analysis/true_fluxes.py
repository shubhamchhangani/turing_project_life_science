"""Authoring tool: the injected ("true") flux map, normalized to glucose uptake = 100.

Net fluxes describe a glycolytic, glutamine-addicted cancer-cell-like phenotype with
appreciable oxidative PPP, malic-enzyme cataplerosis and a small reductive
(reverse IDH) exchange.
"""
NET = {
    "GLCup": 100.0, "PGI": 85.0, "PFK": 93.0, "G6PDH": 15.0,
    "TK1": 4.0, "TA": 4.0, "TK2": 4.0,
    "GAPDH": 190.0, "PK": 186.0, "PHGDH": 4.0,
    "LDH": 167.0, "ALT": 6.0, "PDH": 42.0, "PC": 11.0, "ME": 40.0,
    "CS": 30.0, "IDH": 30.0, "AKGDH": 64.0, "SDH": 64.0, "FH": 64.0, "MDH": 24.0,
    "GOT": 5.0, "GLNup": 40.0, "GLUDH": 34.0,
    "LACout": 167.0, "ALAout": 6.0, "SERout": 4.0, "ASPout": 5.0, "GLUout": 6.0,
    "P5Pout": 3.0, "ACCOAout": 12.0,
}

EXCH = {"PGI": 60.0, "TK1": 15.0, "TA": 25.0, "TK2": 10.0, "IDH": 8.0, "FH": 80.0, "MDH": 40.0, "GLUDH": 50.0}

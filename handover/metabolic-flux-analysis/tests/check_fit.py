"""Predict observations from a submitted map using the independent generator model.

No optimizer, solution import, fitted parameter vector, or reference output is
used. The input is JSON from the verifier's validated CSV parser.
"""
import csv
import json
import math
from pathlib import Path
import sys

import forward_model as model


def weighted_ssr(table):
    data = Path(__file__).with_name("data")
    net = {rid: row["net"] for rid, row in table.items()}
    exchange = {rid: min(row["forward"], row["reverse"]) for rid, row in table.items()}
    predictions = {}
    for experiment in model.EXPERIMENTS:
        dists = model.simulate(net, exchange, model.sources_for(experiment))
        predictions[experiment] = {f: model.raw_mid(dists, f) for f in model.FRAGMENTS}
    residuals = []
    with (data / "mid_measurements.csv").open() as fh:
        for row in csv.DictReader(fh):
            prediction = predictions[row["experiment"]][row["fragment"]][int(row["isotopologue"][1:])]
            residuals.append((prediction - float(row["fraction"])) / float(row["sd"]))
    with (data / "extracellular_rates.csv").open() as fh:
        for row in csv.DictReader(fh):
            if float(row["sd"]) > 0:
                residuals.append((net[row["reaction_id"]] - float(row["rate"])) / float(row["sd"]))
    assert all(math.isfinite(r) for r in residuals), "nonfinite predicted residual"
    return math.fsum(float(r) ** 2 for r in residuals)


if __name__ == "__main__":
    print(weighted_ssr(json.load(sys.stdin)))

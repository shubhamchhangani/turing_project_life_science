"""Verifier for the parallel-tracer 13C-MFA task.

Expected values are the fluxes injected when the synthetic data were generated
(truth by construction), not output of the reference solution. Contract and
internal-consistency checks run first so a failure names the violated field.
"""
import csv
import json
import math
import os
import stat
import re
import subprocess
import sys
from pathlib import Path

import pytest

FLUXES = Path("/root/results/fluxes.csv")
POOLS = Path("/root/results/pools.csv")
NETWORK = Path("/tests/network.tsv")
EXPECTED = json.loads(Path("/tests/expected_fluxes.json").read_text())
DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
HEADER = ["reaction_id", "net_flux", "forward_flux", "reverse_flux"]
BALANCE_TOL = 0.1        # flux units (glucose uptake = 100)
DIRECTION_TOL = 0.05     # |forward - reverse - net|


def network():
    with NETWORK.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [(r["reaction_id"], r["equation"], r["reversible"] == "yes") for r in rows]


def parse_side(side):
    return [m.strip() for m in side.split("+") if m.strip()]


def read_table(path):
    # O_NONBLOCK prevents a FIFO from blocking before its type can be checked.
    fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
    try:
        info = os.fstat(fd)
        assert stat.S_ISREG(info.st_mode), "output must be a regular file"
        assert info.st_size <= 1024 * 1024, "output must not exceed 1 MiB"
        with os.fdopen(fd, encoding="utf-8", newline="") as fh:
            fd = None
            content = fh.read(1024 * 1024 + 1)
        assert len(content) <= 1024 * 1024, "output must not exceed 1 MiB"
        rows = list(csv.reader(content.splitlines()))
        assert rows, "output is empty"
        raw_header = content.split("\n", 1)[0].removesuffix("\r")
        return raw_header, rows[1:]
    finally:
        if fd is not None:
            os.close(fd)


def read_fluxes():
    return read_table(FLUXES)


def pool_table():
    _, rows = read_table(POOLS)
    return {r[0]: float(r[1]) for r in rows}


def test_pool_contract():
    """pools.csv has the exact header, four unknown pools once each, and positive finite decimals."""
    header, rows = read_table(POOLS)
    assert header == "metabolite,pool_size"
    assert len(rows) == 4 and all(len(r) == 2 for r in rows)
    assert sorted(r[0] for r in rows) == sorted(EXPECTED["pools"])
    for met, raw in rows:
        assert DECIMAL.fullmatch(raw), f"{met}: pool size must be an ASCII decimal"
        assert math.isfinite(float(raw)) and float(raw) > 0, f"{met}: pool must be finite and positive"



def flux_table():
    _header, rows = read_fluxes()
    return {r[0]: {"net": float(r[1]), "forward": float(r[2]), "reverse": float(r[3])} for r in rows}


# ----------------------------------------------------------------------------- contract
def test_file_exists():
    """fluxes.csv exists at /root/results/fluxes.csv."""
    assert FLUXES.is_file(), f"missing output file {FLUXES}"


def test_header():
    """Header is exactly reaction_id,net_flux,forward_flux,reverse_flux."""
    header, _ = read_fluxes()
    expected_header = ",".join(HEADER)
    assert header == expected_header, f"header {header!r} != {expected_header!r}"


def test_one_row_per_reaction():
    """Exactly one row for each of the 31 reactions in network.tsv, no extra or duplicate rows."""
    _, rows = read_fluxes()
    ids = [r[0] for r in rows]
    expected = [rid for rid, _, _ in network()]
    assert len(ids) == len(set(ids)), f"duplicate reaction_id rows: {sorted(i for i in ids if ids.count(i) > 1)}"
    assert sorted(ids) == sorted(expected), (
        f"missing {sorted(set(expected) - set(ids))}, unexpected {sorted(set(ids) - set(expected))}")
    assert all(len(r) == 4 for r in rows), "every row must have 4 fields"


def test_values_finite():
    """Flux fields use whitespace-free ASCII decimal notation and have finite values."""
    _, rows = read_fluxes()
    for r in rows:
        for name, raw in zip(HEADER[1:], r[1:]):
            assert DECIMAL.fullmatch(raw), f"{r[0]} {name}: {raw!r} is not an ASCII decimal"
            try:
                v = float(raw)
            except ValueError:
                pytest.fail(f"{r[0]} {name}: '{raw}' is not a number")
            assert math.isfinite(v), f"{r[0]} {name}: {raw} is not finite"


# ----------------------------------------------------------------------------- consistency
def test_glucose_normalization():
    """Glucose uptake (GLCup net flux) is reported as 100."""
    v = flux_table()["GLCup"]["net"]
    assert abs(v - 100.0) <= 1e-6 * 100, f"GLCup net_flux is {v}, must be 100"


def test_directional_consistency():
    """One-way fluxes are nonnegative; direction identity and irreversibility use 0.05 tolerance."""
    f = flux_table()
    for rid, _eq, rev in network():
        row = f[rid]
        assert row["forward"] >= 0 and row["reverse"] >= 0, (
            f"{rid}: directional fluxes must be non-negative ({row})")
        assert abs(row["forward"] - row["reverse"] - row["net"]) <= DIRECTION_TOL, (
            f"{rid}: forward - reverse = {row['forward'] - row['reverse']:.4f} but net = {row['net']:.4f}")
        if not rev:
            assert abs(row["reverse"]) <= DIRECTION_TOL and row["net"] >= -DIRECTION_TOL, (
                f"{rid} is irreversible: reverse_flux must be 0 and net_flux >= 0 ({row})")


def test_mass_balance():
    """Reported net fluxes close the steady-state balance of every intracellular metabolite within 0.1."""
    f = flux_table()
    balance = {}
    for rid, eq, _rev in network():
        lhs, rhs = re.split(r"\s*(?:<=>|->)\s*", eq)
        for m in parse_side(lhs):
            balance[m] = balance.get(m, 0.0) - f[rid]["net"]
        for m in parse_side(rhs):
            balance[m] = balance.get(m, 0.0) + f[rid]["net"]
    bad = {m: round(b, 4) for m, b in balance.items()
           if not m.endswith(".ext") and m != "CO2" and abs(b) > BALANCE_TOL}
    assert not bad, f"unbalanced intracellular metabolites (production - consumption): {bad}"


def test_pool_accuracy():
    """Each inferred pool is within 25 percent of its independently injected value."""
    test_pool_contract()
    for met, value in pool_table().items():
        want = EXPECTED["pools"][met]
        assert abs(value - want) <= 0.25 * want, f"{met}: pool {value}, expected {want} within 25%"


# ----------------------------------------------------------------------------- fluxes
@pytest.mark.parametrize("key", list(EXPECTED["graded"]))
def test_flux_value(key):
    """Graded flux is within its stated relative tolerance of the true (injected) value."""
    spec = EXPECTED["graded"][key]
    got = flux_table()[spec["reaction_id"]][spec["column"]]
    want, rtol = spec["true"], spec["rel_tol"]
    assert abs(got - want) <= rtol * abs(want), (
        f"{key} ({spec['reaction_id']} {spec['column']}): got {got:.3f}, true {want:.3f}, "
        f"allowed relative error {rtol:.0%}")


def test_measurement_fit():
    """The full map predicts all MID fractions and positive-SD rates with weighted SSR <= 2110."""
    # Validate cheap contract checks before starting any scientific computation.
    test_pool_contract()
    test_header()
    test_one_row_per_reaction()
    test_values_finite()
    test_glucose_normalization()
    test_directional_consistency()
    test_mass_balance()
    result = subprocess.run(
        [sys.executable, str(Path(__file__).with_name("check_fit.py"))],
        input=json.dumps({"fluxes": flux_table(), "pools": pool_table()}), capture_output=True, text=True, timeout=60,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )
    assert result.returncode == 0, f"measurement prediction failed: {result.stderr[-2000:]}"
    ssr = float(result.stdout)
    assert math.isfinite(ssr) and ssr <= 2110, f"weighted measurement SSR {ssr:.6f} exceeds 2110"

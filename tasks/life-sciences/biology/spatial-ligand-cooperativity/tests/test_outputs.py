"""Check the complete parameter map and its independent equilibrium predictions."""
import csv
import io
import json
import os
from pathlib import Path
import re
import stat

import numpy as np
import pytest
from physics import IndependentSparseFields

DATA = Path(__file__).parent / 'data'
MODEL = json.loads((DATA / 'model.json').read_text())
TRUTH = json.loads((DATA / 'injected_parameters.json').read_text())
EDGES = [(x, y, d) for y in range(MODEL['height'])
         for x in range(MODEL['width']) for d in ['E', 'N']]
SITES = [p['id'] for p in MODEL['field_sites']]
DECIMAL = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')
INTEGER = re.compile(r'(?:0|[1-9][0-9]*)\Z')


def read_csv(path, header, count):
    descriptor = os.open(path, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW)
    with os.fdopen(descriptor, 'rb') as stream:
        info = os.fstat(stream.fileno())
        assert stat.S_ISREG(info.st_mode), 'Artifact must be a regular file'
        assert info.st_size <= 1048576, 'Artifact exceeds 1 MiB'
        raw = stream.read(1048577)
    assert len(raw) <= 1048576, 'Artifact exceeds 1 MiB'
    text = raw.decode('utf-8')
    assert '\r' not in text.replace('\r\n', ''), 'Use LF or CRLF endings'
    assert text.split('\n', 1)[0].removesuffix('\r') == header, 'Incorrect literal header'
    rows = list(csv.reader(io.StringIO(text, newline=''), strict=True))
    assert len(rows) == count + 1, 'Incorrect row count'
    assert all(len(row) == len(header.split(',')) for row in rows), 'Incorrect column count'
    return rows[1:]


def number(raw, bounds):
    assert DECIMAL.fullmatch(raw), 'Invalid decimal field'
    value = float(raw)
    assert np.isfinite(value), 'Nonfinite parameter'
    assert bounds[0] <= value <= bounds[1], 'Parameter outside stated bounds'
    return value


def contract(directory):
    edges = {}
    for x, y, direction, raw in read_csv(directory / 'couplings.csv', 'x,y,direction,K', len(EDGES)):
        assert INTEGER.fullmatch(x) and INTEGER.fullmatch(y), 'Invalid coordinate syntax'
        key = (int(x), int(y), direction)
        assert key in EDGES and key not in edges, 'Unknown or duplicate edge'
        edges[key] = number(raw, MODEL['coupling_bounds'])
    assert set(edges) == set(EDGES), 'Incomplete edge map'
    fields = {}
    for site, raw in read_csv(directory / 'biases.csv', 'receptor_id,h', len(SITES)):
        assert site in SITES and site not in fields, 'Unknown or duplicate receptor'
        fields[site] = number(raw, MODEL['field_bounds'])
    assert set(fields) == set(SITES), 'Incomplete bias map'
    return np.array([edges[e] for e in EDGES] + [fields[s] for s in SITES])


@pytest.fixture(scope='module')
def submitted():
    return contract(Path('/root/results'))


def test_01_complete_contract(submitted):
    """Both bounded regular CSV files contain every required parameter exactly once."""
    assert len(submitted) == len(EDGES) + len(SITES)


def test_02_every_coupling(submitted):
    """All 646 interactions agree with independently planted generating parameters."""
    assert np.max(np.abs(submitted[:len(EDGES)] - np.array(TRUTH['couplings']))) <= 5e-5


def test_03_every_bias(submitted):
    """All twelve local ligand biases agree with their generating parameters."""
    assert np.max(np.abs(submitted[len(EDGES):] - np.array(TRUTH['biases']))) <= 5e-5


def observed_moments():
    with (DATA / 'edge_moments.csv').open() as stream:
        edges = {(int(r['x']), int(r['y']), r['direction']): float(r['mean_product'])
                 for r in csv.DictReader(stream)}
    with (DATA / 'site_moments.csv').open() as stream:
        fields = {r['receptor_id']: float(r['mean_state']) for r in csv.DictReader(stream)}
    return np.array([edges[e] for e in EDGES] + [fields[s] for s in SITES])


def calculate_moments(parameters):
    sites = [(p['x'], p['y']) for p in MODEL['field_sites']]
    physics = IndependentSparseFields(MODEL['width'], MODEL['height'], sites)
    return physics.observables(parameters)[1]


def test_04_all_equilibrium_moments(submitted):
    """The complete inferred model reproduces every edge and site observation."""
    assert np.max(np.abs(calculate_moments(submitted) - observed_moments())) <= 2e-6

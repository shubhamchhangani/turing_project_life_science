"""Outcome-only, offline verification against a fresh statistical-mechanical calculation."""
import csv
import io
import json
import math
import os
from pathlib import Path
import re
import stat

import numpy as np
import pytest
from physics import DecoratedTorus

RESULTS=Path('/root/results')
DATA=Path(__file__).resolve().parent/'data'
NUMBER=re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')


def read_csv(path,header,count):
    """Bounded, nonblocking read before parsing any submitted values."""
    descriptor=os.open(path,os.O_RDONLY|os.O_NONBLOCK|os.O_NOFOLLOW)
    try:
        info=os.fstat(descriptor)
        assert stat.S_ISREG(info.st_mode), 'Artifact must be a regular file'
        assert 0<info.st_size<=1048576, 'Artifact size outside contract'
        raw=os.read(descriptor,1048577)
        assert len(raw)==info.st_size and len(raw)<=1048576
    finally:
        os.close(descriptor)
    text=raw.decode('utf-8')
    lines=text.splitlines()
    assert lines and lines[0]==','.join(header), 'Incorrect literal CSV header'
    rows=list(csv.reader(io.StringIO(text,newline=''),strict=True))
    assert len(rows)==count+1, 'Incorrect number of rows'
    assert all(len(row)==len(header) for row in rows), 'Incorrect column count'
    return rows[1:]


def number(raw):
    """Validate the numeric grammar before numeric conversion."""
    assert NUMBER.fullmatch(raw), 'Invalid decimal numeric field'
    value=float(raw)
    assert math.isfinite(value), 'Nonfinite numeric field'
    return value


def contract(directory):
    """Validate every output field before any physical comparison."""
    rows=read_csv(directory/'joint.csv',['pattern','probability'],64)
    assert len({r[0] for r in rows})==64, 'Duplicate pattern'
    assert {r[0] for r in rows}=={f'{i:06b}' for i in range(64)}, 'Missing or invalid patterns'
    joint={pattern:number(raw) for pattern,raw in rows}
    assert all(0<=p<=1 for p in joint.values()), 'Probability outside [0,1]'
    summary=read_csv(directory/'thermodynamics.csv',['quantity','value'],1)
    assert summary[0][0]=='log_partition', 'Incorrect thermodynamic quantity'
    logz=number(summary[0][1])
    return np.array([joint[f'{i:06b}'] for i in range(64)]),logz


@pytest.fixture(scope='module')
def submitted():
    return contract(RESULTS)


@pytest.fixture(scope='module')
def calculated(submitted):
    # Depend on the complete contract fixture: malformed files never trigger
    # an expensive calculation or receive credit via another correctness test.
    model=json.loads((DATA/'model.json').read_text())
    physics=DecoratedTorus(model['width'],model['height'])
    with (DATA/'couplings.csv').open() as stream:
        rows=list(csv.DictReader(stream))
    values={(int(r['x']),int(r['y']),0 if r['direction']=='E' else 1):float(r['K']) for r in rows}
    couplings=np.array([values[key] for key in physics.lookup])
    probes=[(p['x'],p['y']) for p in model['probes']]
    return physics.joint(couplings,probes),physics.partition(couplings)


def test_01_complete_contract(submitted):
    """Both regular CSV files must satisfy the complete declared format."""
    assert len(submitted[0])==64


def test_02_normalization(submitted):
    """The unconditional distribution must have total mass one."""
    assert abs(submitted[0].sum()-1)<=2e-7


def test_03_inversion_symmetry(submitted):
    """Zero applied field requires equal probabilities for complementary states."""
    assert np.max(np.abs(submitted[0]-submitted[0][::-1]))<=2e-7


def test_04_all_probe_marginals(submitted):
    """Every probe has activation probability one half in this unbiased model."""
    for position in range(6):
        assert abs(sum(p for i,p in enumerate(submitted[0]) if f'{i:06b}'[position]=='1')-.5)<=2e-7


def test_05_full_joint_distribution(submitted,calculated):
    """Check every pattern, including information in higher-order correlations."""
    assert np.max(np.abs(submitted[0]-calculated[0]))<=1e-7


def test_06_absolute_partition_function(submitted,calculated):
    """Check the absolute finite-system log Z, including all normalization factors."""
    assert abs(submitted[1]-calculated[1])<=1e-7

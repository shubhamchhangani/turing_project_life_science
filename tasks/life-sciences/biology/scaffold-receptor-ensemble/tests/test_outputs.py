"""Validate the complete conditional ensemble against independently generated sums."""
import csv,io,json,math,os,re,stat
from pathlib import Path
import pytest

EXPECTED=json.loads((Path(__file__).parent/'data/expectations.json').read_text())
PATTERNS=[''.join(str((i>>j)&1) for j in range(8)) for i in range(256)]
DECIMAL=re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z',re.ASCII)


def read_rows(path,header,count):
    fd=os.open(path,os.O_RDONLY|os.O_NONBLOCK|os.O_NOFOLLOW)
    with os.fdopen(fd,'rb') as f:
        meta=os.fstat(f.fileno());assert stat.S_ISREG(meta.st_mode),'Output must be a regular file'
        assert meta.st_size<=1048576,'Output exceeds 1 MiB';raw=f.read(1048577)
    assert len(raw)<=1048576,'Output exceeds 1 MiB'
    text=raw.decode('utf-8')
    assert '\r' not in text.replace('\r\n',''),'Invalid line endings'
    assert text.split('\n',1)[0].removesuffix('\r')==header,'Wrong literal header'
    rows=list(csv.reader(io.StringIO(text,newline=''),strict=True))
    assert len(rows)==count+1,'Wrong row count'
    assert all(len(row)==2 for row in rows),'Wrong column count'
    return rows[1:]


def number(raw):
    assert DECIMAL.fullmatch(raw),'Invalid ASCII decimal syntax'
    value=float(raw);assert math.isfinite(value),'Nonfinite value'
    return value


def contract(directory):
    thermo={}
    for key,raw in read_rows(directory/'thermodynamics.csv','metric,value',2):
        assert key in ['log_partition','log_clamp_probability'] and key not in thermo,'Unknown or duplicate metric'
        thermo[key]=number(raw)
    probs={}
    for key,raw in read_rows(directory/'reporter_distribution.csv','pattern,probability',256):
        assert key in PATTERNS and key not in probs,'Unknown or duplicate pattern'
        probs[key]=number(raw);assert 0<=probs[key]<=1,'Probability outside [0,1]'
    assert set(probs)==set(PATTERNS),'Missing pattern'
    return thermo,[probs[p] for p in PATTERNS]


def accuracy(submitted):
    thermo,probs=submitted
    for key,value in thermo.items():assert abs(value-EXPECTED[key])<=2e-6,f'Incorrect {key}'
    assert abs(math.fsum(probs)-1)<=1e-6,'Probabilities do not sum to one'
    for pattern,value,want in zip(PATTERNS,probs,EXPECTED['probabilities']):
        assert abs(value-want)<=5e-8,f'Incorrect joint probability: {pattern}'


@pytest.fixture(scope='module')
def submitted():
    return contract(Path('/root/results'))


def test_01_complete_contract(submitted):
    """Require both bounded regular CSVs and all unique keys with strict decimal syntax."""
    assert len(submitted[0])==2 and len(submitted[1])==256


def test_02_thermodynamics(submitted):
    """Check counting-measure log partition and the joint six-ligand clamp probability."""
    for key,value in submitted[0].items():assert abs(value-EXPECTED[key])<=2e-6,f'Incorrect {key}'


def test_03_normalization(submitted):
    """The full conditional distribution must sum to one within the stated tolerance."""
    assert abs(math.fsum(submitted[1])-1)<=1e-6


def test_04_every_joint_probability(submitted):
    """Check all 256 configurations, preserving spatial ordering and higher correlations."""
    for pattern,value,want in zip(PATTERNS,submitted[1],EXPECTED['probabilities']):
        assert abs(value-want)<=5e-8,f'Incorrect joint probability: {pattern}'

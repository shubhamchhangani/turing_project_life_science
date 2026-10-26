"""Validate all lineage posterior coordinates against independent uniformization."""
import csv,io,json,math,os,re,stat
from pathlib import Path
import pytest

BASE=Path(__file__).parent/'data'
EXPECTED=json.loads((BASE/'expectations.json').read_text())
MODEL=json.loads((BASE/'model.json').read_text())
M=MODEL['rna_capacity'];P=MODEL['protein_capacity'];K=MODEL['count_limit']
TABLES={
    'evidence.csv':('metric,value',['log_evidence'],[EXPECTED['log_evidence']]),
    'division_posterior.csv':('state,probability',[f'{s}:{m}:{p}' for s in range(2) for m in range(M+1) for p in range(P+1)],EXPECTED['parent']),
    'daughter_rna.csv':('counts,probability',[f'{a}:{b}' for a in range(M+1) for b in range(M+1)],[x for row in EXPECTED['daughter_rna'] for x in row]),
    'transcription_counts.csv':('count,probability',[str(k) for k in range(K+1)]+['overflow'],EXPECTED['counts']+[EXPECTED['count_overflow']]),
}
DECIMAL=re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z',re.ASCII)

def contract(directory):
    result={}
    for name,(header,keys,_) in TABLES.items():
        fd=os.open(directory/name,os.O_RDONLY|os.O_NONBLOCK|os.O_NOFOLLOW)
        with os.fdopen(fd,'rb') as f:
            meta=os.fstat(f.fileno());assert stat.S_ISREG(meta.st_mode),'Expected a regular file'
            assert meta.st_size<=1048576,'Output exceeds 1 MiB';raw=f.read(1048577)
        assert len(raw)<=1048576,'Output exceeds 1 MiB'
        text=raw.decode('utf-8');assert '\r' not in text.replace('\r\n',''),'Invalid line endings'
        assert text.split('\n',1)[0].removesuffix('\r')==header,'Wrong literal header'
        rows=list(csv.reader(io.StringIO(text,newline=''),strict=True))
        assert len(rows)==len(keys)+1 and all(len(row)==2 for row in rows),'Wrong CSV dimensions'
        values={}
        for key,raw_value in rows[1:]:
            assert key in keys and key not in values,'Invalid or duplicate key'
            assert DECIMAL.fullmatch(raw_value),'Invalid ASCII decimal'
            v=float(raw_value);assert math.isfinite(v),'Nonfinite value'
            if name!='evidence.csv':assert 0<=v<=1,'Probability outside [0,1]'
            values[key]=v
        assert set(values)==set(keys),'Missing keys'
        result[name]=[values[k] for k in keys]
    return result

def accuracy(data):
    for name,(_,keys,wants) in TABLES.items():
        tol=2e-7 if name=='evidence.csv' else 2e-8
        if name!='evidence.csv':assert abs(math.fsum(data[name])-1)<=1e-6,'Not normalized'
        for k,v,w in zip(keys,data[name],wants):assert abs(v-w)<=tol,f'Incorrect {name} {k}'

@pytest.fixture(scope='module')
def submitted():
    return contract(Path('/root/results'))

def test_01_contract(submitted):
    """Require all four bounded regular CSVs with exact keys and strict numeric fields."""
    assert set(submitted)==set(TABLES)

def test_02_evidence(submitted):
    """Check the full joint observation probability, including integrated photon likelihoods."""
    assert abs(submitted['evidence.csv'][0]-EXPECTED['log_evidence'])<=2e-7

def test_03_normalization(submitted):
    """Check all three complete posterior tables, including the event-count overflow bin."""
    for name in TABLES:
        if name!='evidence.csv':assert abs(math.fsum(submitted[name])-1)<=1e-6

def test_04_every_posterior_coordinate(submitted):
    """Check every mother state, every joint daughter RNA pair and every event count."""
    accuracy(submitted)

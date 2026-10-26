"""Verify the literal artifact, every dosage, and the complete integer MEC cost."""
import csv
import io
import json
import os
from pathlib import Path
import re
import stat
import pytest

DATA=Path(__file__).parent/'data'
MODEL=json.loads((DATA/'model.json').read_text())
with (DATA/'variants.csv').open(newline='') as stream:
    VARIANTS=list(csv.DictReader(stream))
INDEX={row['variant_id']:j for j,row in enumerate(VARIANTS)}
DOSAGES=[int(row['alt_dosage']) for row in VARIANTS]
IDS={f'H{h+1:02d}' for h in range(MODEL['ploidy'])}
ALLELES=re.compile(r'[01]{'+str(MODEL['variants'])+r'}\Z',re.ASCII)


def contract(directory):
    descriptor=os.open(directory/'haplotypes.csv',os.O_RDONLY|os.O_NONBLOCK|os.O_NOFOLLOW)
    with os.fdopen(descriptor,'rb') as stream:
        info=os.fstat(stream.fileno())
        assert stat.S_ISREG(info.st_mode),'Output must be a regular file'
        assert info.st_size<=1048576,'Output exceeds 1 MiB'
        raw=stream.read(1048577)
    assert len(raw)<=1048576,'Output exceeds 1 MiB'
    text=raw.decode('utf-8')
    assert '\r' not in text.replace('\r\n',''),'Use LF or CRLF endings'
    assert text.split('\n',1)[0].removesuffix('\r')=='haplotype_id,alleles','Wrong literal header'
    rows=list(csv.reader(io.StringIO(text,newline=''),strict=True))
    assert len(rows)==MODEL['ploidy']+1,'Wrong number of rows'
    assert all(len(row)==2 for row in rows),'Wrong number of columns'
    found={}
    for name,sequence in rows[1:]:
        assert name in IDS and name not in found,'Unknown or duplicate haplotype ID'
        assert ALLELES.fullmatch(sequence),'Alleles must be exactly 480 ASCII binary digits'
        found[name]=sequence
    assert set(found)==IDS,'Incomplete homolog set'
    return [found[name] for name in sorted(IDS)]


def check_dosages(haplotypes):
    for j,d in enumerate(DOSAGES):
        assert sum(h[j]=='1' for h in haplotypes)==d,f'Incorrect dosage at variant {j+1}'


def correction_cost(haplotypes):
    # Independent pure-Python calculation, with no search code or agent execution.
    by_read={}
    with (DATA/'read_calls.csv').open(newline='') as stream:
        for row in csv.DictReader(stream):
            by_read.setdefault(row['read_id'],[]).append((INDEX[row['variant_id']],row['allele'],int(row['weight'])))
    assert len(by_read)==MODEL['reads'],'Verifier data integrity error'
    total=0
    for calls in by_read.values():
        total+=min(sum(w for j,a,w in calls if h[j]!=a) for h in haplotypes)
    return total


@pytest.fixture(scope='module')
def submitted():
    return contract(Path('/root/results'))


def test_01_complete_contract(submitted):
    """All sixteen homolog sequences satisfy the full literal output contract."""
    assert len(submitted)==MODEL['ploidy']


def test_02_all_variant_dosages(submitted):
    """Every variant has exactly its supplied alternate-allele copy number."""
    check_dosages(submitted)


def test_03_all_read_corrections(submitted):
    """Every observed read contributes its complete best-homolog weighted distance."""
    check_dosages(submitted)
    assert correction_cost(submitted)<=MODEL['maximum_weighted_corrections'],'Correction budget exceeded'

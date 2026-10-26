"""Check all Bayesian quantities against independent partition-ratio expectations."""
import csv,io,json,math,os,re,stat
from pathlib import Path
import pytest
DATA=Path(__file__).parent/'data'
MODEL=json.loads((DATA/'model.json').read_text())
EXPECTED=json.loads((DATA/'expectations.json').read_text())
READS=[r['id'] for r in MODEL['reads']]
assert len(EXPECTED['read_expectations'])==len(READS)
assert all(len(row)==2 for row in EXPECTED['read_expectations'])
DECIMAL=re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z',re.ASCII)


def read_rows(path,header,count):
 fd=os.open(path,os.O_RDONLY|os.O_NONBLOCK|os.O_NOFOLLOW)
 with os.fdopen(fd,'rb') as f:
  meta=os.fstat(f.fileno());assert stat.S_ISREG(meta.st_mode),'Output must be a regular file'
  assert meta.st_size<=1048576,'Output exceeds 1 MiB';raw=f.read(1048577)
 assert len(raw)<=1048576,'Output exceeds 1 MiB';text=raw.decode('utf-8')
 assert '\r' not in text.replace('\r\n',''),'Use LF or CRLF line endings'
 assert text.split('\n',1)[0].removesuffix('\r')==header,'Wrong literal header'
 rows=list(csv.reader(io.StringIO(text,newline=''),strict=True))
 assert len(rows)==count+1,'Wrong row count'
 assert all(len(row)==len(header.split(',')) for row in rows),'Wrong column count'
 return rows[1:]


def number(raw,lo=None,hi=None):
 assert DECIMAL.fullmatch(raw),'Invalid ASCII decimal syntax'
 x=float(raw);assert math.isfinite(x),'Nonfinite value'
 assert lo is None or x>=lo,'Below lower bound'
 assert hi is None or x<=hi,'Above upper bound'
 return x


def contract(directory):
 rows=read_rows(directory/'evidence.csv','metric,value',1)
 assert rows[0][0]=='log_evidence','Unknown evidence metric'
 z=number(rows[0][1]);found={}
 for name,p,m in read_rows(directory/'read_posteriors.csv','read_id,p_error_free,expected_errors',len(READS)):
  assert name in READS and name not in found,'Unknown or duplicate read ID'
  found[name]=(number(p,0,1),number(m,0,8))
 assert set(found)==set(READS),'Missing read'
 return z,[found[r] for r in READS]


def accuracy(submitted):
 z,rows=submitted
 assert abs(z-EXPECTED['log_evidence'])<=2e-6,'Incorrect log evidence'
 for actual,want in zip(rows,EXPECTED['read_expectations']):
  assert all(abs(a-b)<=5e-7 for a,b in zip(actual,want)),'Incorrect read posterior'


@pytest.fixture(scope='module')
def submitted():
 return contract(Path('/root/results'))


def test_01_complete_contract(submitted):
 """Both regular bounded CSVs contain all finite, correctly formatted quantities."""
 assert len(submitted[1])==114


def test_02_normalized_log_evidence(submitted):
 """The evidence includes all phase and origin prior normalizers."""
 assert abs(submitted[0]-EXPECTED['log_evidence'])<=2e-6


def test_03_every_read_posterior(submitted):
 """Both posterior quantities for every read condition on the complete dataset."""
 for read,actual,want in zip(READS,submitted[1],EXPECTED['read_expectations']):
  assert abs(actual[0]-want[0])<=5e-7,f'Incorrect error-free probability: {read}'
  assert abs(actual[1]-want[1])<=5e-7,f'Incorrect expected error count: {read}'

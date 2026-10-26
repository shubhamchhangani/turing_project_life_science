"""Check declared artifacts against independent injected targets and equations."""
import csv,json,math,os,re,stat,subprocess,sys
from pathlib import Path
import pytest

EXPECTED=json.loads(Path('/tests/expected_fluxes.json').read_text())
DECIMAL=re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?')
POPS=('A','B')
RESULTS=Path('/root/results')
def network():
    with Path('/tests/network.tsv').open() as f:return list(csv.DictReader(f,delimiter='\t'))

def read_table(name):
    fd=os.open(RESULTS/name,os.O_RDONLY|os.O_NONBLOCK)
    try:
        info=os.fstat(fd)
        assert stat.S_ISREG(info.st_mode),'output must be a regular file'
        assert info.st_size<=1024*1024,'output exceeds 1 MiB'
        with os.fdopen(fd,encoding='utf-8',newline='') as f:
            fd=None;content=f.read(1024*1024+1)
        assert len(content)<=1024*1024,'output exceeds 1 MiB'
        rows=list(csv.reader(content.splitlines()))
        assert rows,'empty output'
        return content.split('\n',1)[0].removesuffix('\r'),rows[1:]
    finally:
        if fd is not None:os.close(fd)

def finite_decimal(raw):
    assert DECIMAL.fullmatch(raw),f'not an ASCII decimal: {raw!r}'
    value=float(raw)
    assert math.isfinite(value),f'non-finite number: {raw!r}'
    return value

def flux_table():
    _,rows=read_table('fluxes.csv')
    out={p:{} for p in POPS}
    for pop,rid,n,f,b in rows:out[pop][rid]={'net':float(n),'forward':float(f),'reverse':float(b)}
    return out

def pool_table():
    _,rows=read_table('pools.csv')
    out={p:{} for p in POPS}
    for pop,met,value in rows:out[pop][met]=float(value)
    return out

def test_flux_contract():
    """Flux file has a literal header, exactly 62 unique population/reaction rows, and finite ASCII numbers."""
    header,rows=read_table('fluxes.csv')
    assert header=='population,reaction_id,net_flux,forward_flux,reverse_flux'
    assert len(rows)==62 and all(len(r)==5 for r in rows)
    want=sorted((p,r['reaction_id']) for p in POPS for r in network())
    assert sorted((r[0],r[1]) for r in rows)==want
    for row in rows:
        for raw in row[2:]:finite_decimal(raw)

def test_pool_contract():
    """Pool file has a literal header, exactly eight unique required pairs, and positive finite ASCII numbers."""
    header,rows=read_table('pools.csv')
    assert header=='population,metabolite,pool_size'
    assert len(rows)==8 and all(len(r)==3 for r in rows)
    assert sorted((r[0],r[1]) for r in rows)==sorted((p,m) for p in POPS for m in EXPECTED['pools'][p])
    for row in rows:assert finite_decimal(row[2])>0,'pool size must be positive'

def test_glucose_normalization():
    """Each population has glucose uptake 100 within absolute tolerance 0.0001."""
    values=flux_table()
    for pop in POPS:assert abs(values[pop]['GLCup']['net']-100)<=.0001

def test_directional_consistency():
    """Every one-way rate is nonnegative; direction identity and irreversibility use 0.05 tolerance."""
    values=flux_table()
    for pop in POPS:
        for reaction in network():
            r=values[pop][reaction['reaction_id']]
            assert r['forward']>=0 and r['reverse']>=0
            assert abs(r['forward']-r['reverse']-r['net'])<=.05
            if reaction['reversible']!='yes':assert r['reverse']<=.05 and r['net']>=-.05

def test_mass_balance():
    """Both populations separately close every intracellular steady-state balance within 0.1."""
    values=flux_table()
    for pop in POPS:
        balance={}
        for reaction in network():
            left,right=re.split(r'\s*(?:<=>|->)\s*',reaction['equation'])
            for sign,side in [(-1,left),(1,right)]:
                for item in side.split('+'):
                    met=item.strip();balance[met]=balance.get(met,0)+sign*values[pop][reaction['reaction_id']]['net']
        bad={m:v for m,v in balance.items() if not m.endswith('.ext') and m!='CO2' and abs(v)>.1}
        assert not bad,f'{pop}: unbalanced metabolites {bad}'

@pytest.mark.parametrize('key',list(EXPECTED['graded']))
def test_flux_accuracy(key):
    """Each of six flux targets per population matches independently injected truth within the stated tolerance."""
    spec=EXPECTED['graded'][key]
    got=flux_table()[spec['population']][spec['reaction_id']][spec['column']]
    assert abs(got-spec['true'])<=spec['rel_tol']*abs(spec['true']),f'{key}: {got} outside declared tolerance'

def test_pool_accuracy():
    """Every unknown pool matches its population-specific injected amount within 25 percent."""
    test_pool_contract();values=pool_table()
    for pop in POPS:
        for met,want in EXPECTED['pools'][pop].items():
            assert abs(values[pop][met]-want)<=.25*want,f'{pop}/{met}: incorrect pool size'

def test_measurement_fit():
    """All four tracer programs and positive-SD bulk rates have total weighted SSR at most 4090."""
    test_flux_contract();test_pool_contract();test_glucose_normalization()
    test_directional_consistency();test_mass_balance()
    result=subprocess.run([sys.executable,'/tests/check_fit.py'],
        input=json.dumps({'fluxes':flux_table(),'pools':pool_table()}),capture_output=True,text=True,timeout=60,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    assert result.returncode==0,f'prediction failed: {result.stderr[-1500:]}'
    ssr=float(result.stdout)
    assert math.isfinite(ssr) and ssr<=4090,f'weighted SSR {ssr:.6f} exceeds 4090'

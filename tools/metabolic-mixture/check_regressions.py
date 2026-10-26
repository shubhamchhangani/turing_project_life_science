"""Author-only adversarial artifacts for the population-deconvolution verifier."""
import csv,json,os,subprocess
from pathlib import Path
import sys
sys.path.insert(0,'/tests')
from network_model import FluxModel

truth=json.loads(Path('/author/injected_truth.json').read_text())
with Path('/tests/data/network.tsv').open() as f:network=list(csv.DictReader(f,delimiter='\t'))
fm=FluxModel([r['reaction_id'] for r in network],{r['reaction_id']:r['reversible']=='yes' for r in network})
def flux_rows():
    rows=[]
    for pop,spec in truth.items():
        for rid in fm.rxns:
            v=spec['net'][rid];x=spec['exchange'].get(rid,0)
            rows.append([pop,rid,v,max(v,0)+x,max(-v,0)+x])
    return rows
def pool_rows():return [[p,m,s['pools'][m]] for p,s in truth.items() for m in ['OAA','Cit','AKG','Mal']]
headers={'fluxes.csv':['population','reaction_id','net_flux','forward_flux','reverse_flux'],'pools.csv':['population','metabolite','pool_size']}
root=Path('/root/results');root.mkdir(exist_ok=True)
results=[]
def check(label,mutate=None,file='fluxes.csv',special=None,transform=None,valid=False):
    for name,rows in [('fluxes.csv',flux_rows()),('pools.csv',pool_rows())]:
        p=root/name;p.unlink(missing_ok=True)
        h=headers[name][:]
        if name==file and mutate:mutate(h,rows)
        with p.open('w',newline='') as f:
            w=csv.writer(f);w.writerow(h);w.writerows(rows)
    p=root/file
    if special:
        p.unlink()
        if special=='fifo':os.mkfifo(p)
        elif special=='oversize':p.write_bytes(b'x'*(1024*1024+1))
        elif special=='empty':p.write_bytes(b'')
    if transform:p.write_bytes(transform(p.read_bytes()))
    test=subprocess.run(['pytest','-q','/tests/test_outputs.py'],capture_output=True,text=True,timeout=20)
    passed=test.returncode==0
    assert passed==valid,(label,test.stdout,test.stderr)
    results.append({'case':label,'expected_valid':valid,'actual_exit_code':test.returncode})
    print(label,'ACCEPTED' if passed else 'REJECTED',flush=True)

check('injected truth',valid=True)
check('literal LF',transform=lambda b:b.replace(b'\r\n',b'\n'),valid=True)
check('ASCII signed exponent',mutate=lambda h,r:r[0].__setitem__(2,'+1e2'),valid=True)
check('row ordering',mutate=lambda h,r:r.reverse(),valid=True)
check('header whitespace',mutate=lambda h,r:h.__setitem__(0,'population '))
check('quoted literal header',transform=lambda b:b.replace(b'population,',b'"population",',1))
check('missing row',mutate=lambda h,r:r.pop())
check('duplicate row',mutate=lambda h,r:r.__setitem__(-1,r[0][:]))
check('unknown population',mutate=lambda h,r:r[0].__setitem__(0,'C'))
check('leading numeric space',mutate=lambda h,r:r[0].__setitem__(2,' 100'))
check('trailing numeric space',mutate=lambda h,r:r[0].__setitem__(2,'100 '))
check('Unicode digits',mutate=lambda h,r:r[0].__setitem__(2,'１００'))
check('underscore',mutate=lambda h,r:r[0].__setitem__(2,'1_00'))
check('NaN',mutate=lambda h,r:r[0].__setitem__(2,'nan'))
check('infinite decimal',mutate=lambda h,r:r[0].__setitem__(2,'1e309'))
check('negative directional rate',mutate=lambda h,r:r[0].__setitem__(slice(3,5),[99.999,-.001]))
def inflate_pgi(h,rows):
    r=next(r for r in rows if r[:2]==['A','PGI']);r[3]+=1000;r[4]+=1000
check('PGI exchange inflation preserves all net fluxes and graded targets',mutate=inflate_pgi)
def compensate_ppp(h,rows):
    for pop,shift in [('A',6),('B',-6*.37/.63)]:
        z=[truth[pop]['net'][r] for r in fm.free];z[fm.free.index('G6PDH')]+=shift
        net=fm.net(z)
        for row in rows:
            if row[0]==pop:
                v=net[row[1]];x=truth[pop]['exchange'].get(row[1],0)
                row[2:]=[v,max(v,0)+x,max(-v,0)+x]
check('opposite PPP shifts preserve balance bulk rates and truth tolerances',mutate=compensate_ppp)
for special in ('missing','empty','oversize','fifo'):check('flux '+special,special=special)
check('pool header whitespace',file='pools.csv',mutate=lambda h,r:h.__setitem__(2,'pool_size '))
check('pool duplicate',file='pools.csv',mutate=lambda h,r:r.__setitem__(-1,r[0][:]))
check('pool unknown population',file='pools.csv',mutate=lambda h,r:r[0].__setitem__(0,'C'))
check('pool numeric whitespace',file='pools.csv',mutate=lambda h,r:r[0].__setitem__(2,' 8'))
check('pool zero',file='pools.csv',mutate=lambda h,r:r[0].__setitem__(2,'0'))
check('pool wrong amount',file='pools.csv',mutate=lambda h,r:r[0].__setitem__(2,'1'))
for special in ('missing','empty','oversize','fifo'):check('pool '+special,file='pools.csv',special=special)
Path('/author_reports/regression-results.json').write_text(json.dumps(results,indent=2)+'\n')

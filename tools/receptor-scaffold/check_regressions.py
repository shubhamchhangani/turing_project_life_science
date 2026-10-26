"""Exercise format boundaries and every independently graded output coordinate."""
import csv,importlib.util,json,os,shutil,tempfile
from pathlib import Path

root=Path('tasks/life-sciences/biology/scaffold-receptor-ensemble')
spec=importlib.util.spec_from_file_location('verifier',root/'tests/test_outputs.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
expected=v.EXPECTED;results=[]
with tempfile.TemporaryDirectory() as tmp:
    out=Path(tmp)
    def baseline():
        for p in out.iterdir():p.unlink()
        with (out/'thermodynamics.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['metric','value'])
            for key in ['log_partition','log_clamp_probability']:w.writerow([key,format(expected[key],'.16g')])
        with (out/'reporter_distribution.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['pattern','probability'])
            w.writerows((p,format(x,'.16g')) for p,x in zip(v.PATTERNS,expected['probabilities']))
    def check(name,valid):
        try:v.accuracy(v.contract(out));accepted=True
        except (AssertionError,OSError,ValueError,UnicodeError,csv.Error):accepted=False
        assert accepted==valid,(name,accepted)
        results.append({'case':name,'expected_accept':valid})
    baseline();check('reference',True)
    for raw in [' 0.1','0.1 ','NaN','inf','1_0','０.１','1e999','-0.1','1.1','', '0x1p-2']:
        baseline();p=out/'reporter_distribution.csv';rows=p.read_text().splitlines();rows[1]='00000000,'+raw;p.write_text('\n'.join(rows)+'\n');check('numeric '+repr(raw),False)
    for i in range(256):
        baseline();p=out/'reporter_distribution.csv';rows=p.read_text().splitlines()
        for j,delta in [(i,1e-5),((i+1)%256,-1e-5)]:rows[j+1]=v.PATTERNS[j]+','+str(expected['probabilities'][j]+delta)
        p.write_text('\n'.join(rows)+'\n');check(f'wrong joint probability {i} with preserved normalization',False)
    for key in ['log_partition','log_clamp_probability']:
        baseline();p=out/'thermodynamics.csv';s=p.read_text();s=s.replace(key+','+format(expected[key],'.16g'),key+','+str(expected[key]+1e-3));p.write_text(s);check('wrong '+key,False)
    baseline()
    marginal=[sum(p for i,p in enumerate(expected['probabilities']) if i>>j&1) for j in range(8)]
    independent=[]
    for i in range(256):
        product=1.0
        for j,p in enumerate(marginal):product*=p if i>>j&1 else 1-p
        independent.append(product)
    with (out/'reporter_distribution.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['pattern','probability']);w.writerows(zip(v.PATTERNS,independent))
    check('product law with correct one-site marginals and normalization',False)
    for mode in ['duplicate','missing','extra','extra_column','header','symlink','fifo','directory','oversize','bad_utf8']:
        baseline();p=out/'reporter_distribution.csv';text=p.read_text();rows=text.splitlines()
        if mode=='duplicate':rows[-1]=rows[1];p.write_text('\n'.join(rows)+'\n')
        elif mode=='missing':p.write_text('\n'.join(rows[:-1])+'\n')
        elif mode=='extra':p.write_text(text+rows[1]+'\n')
        elif mode=='extra_column':p.write_text(text.replace(rows[1],rows[1]+',x'))
        elif mode=='header':p.write_text(text.replace('pattern,probability','pattern, probability'))
        elif mode=='symlink':p.unlink();p.symlink_to(out/'thermodynamics.csv')
        elif mode=='fifo':p.unlink();os.mkfifo(p)
        elif mode=='directory':p.unlink();p.mkdir()
        elif mode=='oversize':p.write_bytes(b'x'*1048577)
        elif mode=='bad_utf8':p.write_bytes(b'\xff')
        check(mode,False)
        if p.is_dir():p.rmdir()
    baseline()
    for p in out.iterdir():
        lines=p.read_text().splitlines();p.write_bytes(('\r\n'.join([lines[0]]+list(reversed(lines[1:])))+'\r\n').encode())
    check('reordered CRLF',True)
    baseline()
    for p in out.iterdir():
        lines=p.read_text().splitlines();p.write_text('\n'.join([lines[0]]+[a+','+format(float(b),'.15e') for a,b in (row.split(',') for row in lines[1:])])+'\n')
    check('scientific notation',True)
report={'passed':True,'valid_cases':sum(r['expected_accept'] for r in results),
        'invalid_cases':sum(not r['expected_accept'] for r in results),'cases':results}
Path('tools/receptor-scaffold/draft/regressions.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k!='cases'})

"""Exercise every numerical coordinate and format boundary before model trials."""
import csv,importlib.util,json,os,tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];AUTHOR=Path(__file__).parent
task=ROOT/'tasks/life-sciences/biology/lineage-photon-inference'
spec=importlib.util.spec_from_file_location('verifier',task/'tests/test_outputs.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
cases=[]

with tempfile.TemporaryDirectory() as tmp:
    out=Path(tmp)
    def baseline():
        for p in out.iterdir():
            if p.is_dir() and not p.is_symlink():p.rmdir()
            else:p.unlink()
        for name,(header,keys,values) in v.TABLES.items():
            (out/name).write_text(header+'\n'+''.join(f'{k},{x:.16g}\n' for k,x in zip(keys,values)))
    def check(name,valid):
        try:v.accuracy(v.contract(out));accepted=True
        except (AssertionError,ValueError,OSError,UnicodeError,csv.Error):accepted=False
        assert accepted==valid,(name,accepted)
        cases.append({'case':name,'expected_accept':valid})
    baseline();check('independent expectations',True)
    for name,(_,keys,values) in v.TABLES.items():
        for i,key in enumerate(keys):
            baseline();p=out/name;rows=p.read_text().splitlines();delta=1e-4
            if name=='evidence.csv':rows[i+1]=key+','+str(values[i]+delta)
            else:
                # Move positive mass from a sufficiently large bin, preserving normalization.
                j=max((j for j in range(len(values)) if j!=i),key=lambda j:values[j])
                assert values[j]>delta
                rows[i+1]=key+','+str(values[i]+delta);rows[j+1]=keys[j]+','+str(values[j]-delta)
            p.write_text('\n'.join(rows)+'\n');check(f'wrong {name} {key}',False)
        for raw in [' 0.1','0.1 ','NaN','inf','1_0','０.１','1e999','','0x1p-2']:
            baseline();p=out/name;rows=p.read_text().splitlines();rows[1]=keys[0]+','+raw;p.write_text('\n'.join(rows)+'\n');check(f'{name} decimal {raw!r}',False)
        for mode in ['duplicate','missing','extra','column','header','symlink','fifo','directory','oversize','utf8']:
            if mode=='duplicate' and len(keys)==1:continue
            baseline();p=out/name;s=p.read_text();rows=s.splitlines()
            if mode=='duplicate':rows[-1]=rows[1];p.write_text('\n'.join(rows)+'\n')
            elif mode=='missing':p.write_text('\n'.join(rows[:-1])+'\n')
            elif mode=='extra':p.write_text(s+rows[1]+'\n')
            elif mode=='column':p.write_text(s.replace(rows[1],rows[1]+',x'))
            elif mode=='header':p.write_text('wrong,header\n'+'\n'.join(rows[1:])+'\n')
            elif mode=='symlink':p.unlink();p.symlink_to(out/'unavailable.csv')
            elif mode=='fifo':p.unlink();os.mkfifo(p)
            elif mode=='directory':p.unlink();p.mkdir()
            elif mode=='oversize':p.write_bytes(b'x'*1048577)
            elif mode=='utf8':p.write_bytes(b'\xff')
            check(f'{name} {mode}',False)
    baseline()
    for p in out.iterdir():
        lines=p.read_text().splitlines();p.write_bytes(('\r\n'.join([lines[0]]+lines[:0:-1])+'\r\n').encode())
    check('reordered CRLF',True)
    baseline()
    for p in out.iterdir():
        lines=p.read_text().splitlines();p.write_text(lines[0]+'\n'+''.join(f'{k},{float(x):.15e}\n' for k,x in (row.split(',') for row in lines[1:])))
    check('scientific notation',True)
    # Correct marginals do not substitute for the conserved joint daughter law.
    baseline();name='daughter_rna.csv';header,keys,vals=v.TABLES[name];n=v.M+1
    aa=[sum(vals[i*n:(i+1)*n]) for i in range(n)];bb=[sum(vals[j::n]) for j in range(n)]
    (out/name).write_text(header+'\n'+''.join(f'{i}:{j},{aa[i]*bb[j]:.16g}\n' for i in range(n) for j in range(n)))
    check('independent daughters with correct marginals',False)
    baseline();name='transcription_counts.csv';header,keys,vals=v.TABLES[name];z=sum(vals[:-1])
    (out/name).write_text(header+'\n'+''.join(f'{k},{x/z:.16g}\n' for k,x in zip(keys[:-1],vals[:-1]))+'overflow,0\n')
    check('conditioning away event-count overflow',False)

report=dict(passed=True,valid_cases=sum(c['expected_accept'] for c in cases),invalid_cases=sum(not c['expected_accept'] for c in cases),cases=cases)
(AUTHOR/'draft/regressions.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k!='cases'})

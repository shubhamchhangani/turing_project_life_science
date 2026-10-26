"""Run inside the verifier image with author tools mounted at /author.

These are author-side regression cases, excluded from the task package.
"""
import csv, io, os, pathlib, runpy, subprocess
truth = runpy.run_path('/author/true_fluxes.py')
header = ['reaction_id', 'net_flux', 'forward_flux', 'reverse_flux']
rows = [[rid, v, max(v, 0) + truth['EXCH'].get(rid, 0), max(-v, 0) + truth['EXCH'].get(rid, 0)] for rid,v in truth['NET'].items()]
p = pathlib.Path('/root/results/fluxes.csv')
p.parent.mkdir(parents=True, exist_ok=True)
def check(label, change=None, special=None, valid=False, transform=None):
    p.unlink(missing_ok=True)
    out = [r[:] for r in rows]
    h = header[:]
    if change: change(h,out)
    if special == 'fifo': os.mkfifo(p)
    elif special == 'empty': p.write_text('')
    elif special == 'oversize': p.write_bytes(b'x' * (1024*1024+1))
    elif special != 'missing':
        with p.open('w', newline='') as f:
            w=csv.writer(f); w.writerow(h); w.writerows(out)
    if transform:
        p.write_bytes(transform(p.read_bytes()))
    result = subprocess.run(['pytest','-q','/tests/test_outputs.py'],capture_output=True,text=True,timeout=75)
    assert (result.returncode == 0) == valid, (label,result.stdout,result.stderr)
    print(label, 'PASS' if valid else 'REJECTED', flush=True)
check('injected truth',valid=True)
check('header whitespace',lambda h,r: h.__setitem__(0, 'reaction_id '))
check('negative directional flux',lambda h,r: r[0].__setitem__(slice(2,4), [r[0][1]-0.001,-0.001]))
check('missing row',lambda h,r: r.pop())
check('duplicate row',lambda h,r: r.append(r[0][:]))
check('nonfinite value',lambda h,r: r[0].__setitem__(1,'nan'))
check('inaccurate flux',lambda h,r: next(x for x in r if x[0]=='CS').__setitem__(1,100))
check('FIFO',special='fifo')
check('empty file',special='empty')
check('oversize file',special='oversize')
check('missing file',special='missing')

check('literal header LF',valid=True,transform=lambda b: b.replace(b'\r\n', b'\n'))
check('literal header CRLF',valid=True)
check('quoted first header field',transform=lambda b: b.replace(b'reaction_id,', b'"reaction_id",', 1))
check('fully quoted header',transform=lambda b: b.replace(b'reaction_id,net_flux,forward_flux,reverse_flux', b'"reaction_id","net_flux","forward_flux","reverse_flux"', 1))
check('header leading whitespace',transform=lambda b: b' '+b)
check('header trailing whitespace',transform=lambda b: b.replace(b'reverse_flux\r\n', b'reverse_flux \r\n', 1))

check('leading numeric space', lambda h,r: r[0].__setitem__(1,' 100'))
check('trailing numeric space', lambda h,r: r[0].__setitem__(1,'100 '))
check('numeric underscore', lambda h,r: r[0].__setitem__(1,'1_00'))
check('unicode digits', lambda h,r: r[0].__setitem__(1,'１００'))
check('decimal exponent', lambda h,r: r[0].__setitem__(1,'1e2'),valid=True)
check('signed decimal', lambda h,r: r[0].__setitem__(1,'+100.0'),valid=True)
check('overflowing decimal', lambda h,r: r[0].__setitem__(1,'1e999'))
def alter_exchange(r, rid, amount):
    row = next(x for x in r if x[0] == rid)
    row[2] += amount
    row[3] += amount
check('TK1 exchange corruption',lambda h,r: alter_exchange(r,'TK1',1000))
check('TK2 exchange corruption',lambda h,r: alter_exchange(r,'TK2',1000))
check('GLUDH exchange corruption',lambda h,r: alter_exchange(r,'GLUDH',1000))

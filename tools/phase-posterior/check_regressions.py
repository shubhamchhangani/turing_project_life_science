"""Exercise the complete output contract and all posterior accuracy requirements."""
import csv,importlib.util,io,json,math,os,tempfile
from pathlib import Path
ROOT=Path('tasks/life-sciences/biology/polyploid-phase-posterior')
spec=importlib.util.spec_from_file_location('posterior_verifier',ROOT/'tests/test_outputs.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
z=v.EXPECTED['log_evidence'];a=v.EXPECTED['read_expectations']
rows=[[name,f'{x:.12f}',f'{m:.12f}'] for name,(x,m) in zip(v.READS,a)]

def render(header,rows,ending='\n',quoting=csv.QUOTE_MINIMAL):
 f=io.StringIO(newline='');f.write(header+ending);w=csv.writer(f,lineterminator=ending,quoting=quoting);w.writerows(rows);return f.getvalue().encode()

def representation(rs=None,zz=None,ending='\n',quoting=csv.QUOTE_MINIMAL):
 return {'evidence.csv':render('metric,value',[['log_evidence',f'{z:.12f}' if zz is None else zz]],ending,quoting),'read_posteriors.csv':render('read_id,p_error_free,expected_errors',rows if rs is None else rs,ending,quoting)}

def changed(i,j,value):
 copy=[r[:] for r in rows];copy[i][j]=value;return representation(copy)

cases=[]
with tempfile.TemporaryDirectory(prefix='posterior-regression-') as tmp:
 d=Path(tmp)
 def attempt(name,files,accepted):
  for p in d.iterdir():
   if p.is_dir() and not p.is_symlink():p.rmdir()
   else:p.unlink()
  for name2,data in files.items():(d/name2).write_bytes(data)
  try:v.accuracy(v.contract(d));ok=True
  except (AssertionError,OSError,ValueError,csv.Error):ok=False
  assert ok==accepted,(name,ok,accepted)
  cases.append(dict(name=name,expected='accept' if accepted else 'reject',passed=True))
 valid={'decimal':representation(),'CRLF':representation(ending='\r\n'),'quoted_fields':representation(quoting=csv.QUOTE_ALL),'reversed_reads':representation(rows[::-1]),'scientific':representation([[r,f'{x:+.12e}',f'{m:+.12E}'] for r,(x,m) in zip(v.READS,a)],f'{z:+.12e}'),'nine_decimals':representation([[r,f'{x:.9f}',f'{m:.9f}'] for r,(x,m) in zip(v.READS,a)],f'{z:.9f}')}
 for name,files in valid.items():attempt(name,files,True)
 for bad in [' '+rows[0][1],rows[0][1]+' ','NaN','inf','-inf','1e999','0_5','０.５','.','+','0x1p0','1,000','-0.0001','1.0001']:
  attempt('bad_probability_'+bad,changed(0,1,bad),False)
 for bad in [' 0','NaN','1e999','9','-1']:
  attempt('bad_expected_errors_'+bad,changed(0,2,bad),False)
 for bad in [' '+str(z),'NaN','1e999',f'{z+114*math.log(6):.10f}']:
  attempt('bad_evidence_'+bad,representation(zz=bad),False)
 for i in [0,56,113]:
  for column in [1,2]:attempt(f'accuracy_read_{i}_column_{column}',changed(i,column,f'{float(rows[i][column])+1e-4:.12f}'),False)
 attempt('wrong_evidence',representation(zz=f'{z+.001:.10f}'),False)
 for name,rs in [('missing_read',rows[:-1]),('duplicate_read',rows[:-1]+[rows[0]]),('extra_read',rows+[rows[0]]),('extra_column',[rows[0]+['x']]+rows[1:])]:attempt(name,representation(rs),False)
 attempt('unknown_read',changed(0,0,'r9999'),False)
 attempt('whitespace_read_id',changed(0,0,' r0001'),False)
 baseline=[[r['id'],str(math.prod(1-c['error_probability'] for c in r['calls'])),str(sum(c['error_probability'] for c in r['calls']))] for r in v.MODEL['reads']]
 attempt('independent_read_posterior',representation(baseline),False)
 for filename in ['evidence.csv','read_posteriors.csv']:
  for name,data in [('empty',b''),('too_large',b'x'*1048577),('invalid_utf8',b'\xff'),('BOM',b'\xef\xbb\xbf'+representation()[filename]),('blank_row',representation()[filename]+b'\n'),('bare_CR',representation()[filename].replace(b'\n',b'\r'))]:
   fs=representation();fs[filename]=data;attempt(filename+'_'+name,fs,False)
  for kind in ['missing','symlink','fifo','directory']:
   fs=representation();del fs[filename]
   for p in d.iterdir():p.unlink()
   for name,data in fs.items():(d/name).write_bytes(data)
   p=d/filename
   if kind=='symlink':(d/'other').write_bytes(representation()[filename]);p.symlink_to(d/'other')
   elif kind=='fifo':os.mkfifo(p)
   elif kind=='directory':p.mkdir()
   try:v.accuracy(v.contract(d))
   except (AssertionError,OSError,ValueError,csv.Error):cases.append(dict(name=filename+'_'+kind,expected='reject',passed=True))
   else:raise AssertionError(kind)
   if kind=='directory':p.rmdir()
   elif kind!='missing':p.unlink()
Path('tools/phase-posterior/draft/regressions.json').write_text(json.dumps(dict(passed=True,cases=cases),indent=2)+'\n')
print(f'{len(valid)} valid and {len(cases)-len(valid)} invalid cases passed.')

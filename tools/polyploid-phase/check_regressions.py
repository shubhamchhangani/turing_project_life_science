"""Exercise valid alternative representations and adversarial artifacts end to end."""
import csv,importlib.util,io,json,os,tempfile
from pathlib import Path
ROOT=Path('tasks/life-sciences/biology/polyploid-haplotype-assembly')
spec=importlib.util.spec_from_file_location('phase_verifier',ROOT/'tests/test_outputs.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
truth=json.loads(Path('tools/polyploid-phase/draft/generation.json').read_text())['original_haplotypes']
rows=[(f'H{i+1:02d}',h) for i,h in enumerate(truth)]

def render(rows,ending='\n',quoting=csv.QUOTE_MINIMAL):
 f=io.StringIO(newline='');f.write('haplotype_id,alleles'+ending);w=csv.writer(f,lineterminator=ending,quoting=quoting);w.writerows(rows);return f.getvalue().encode()

def full(directory):
 h=v.contract(directory);v.check_dosages(h);assert v.correction_cost(h)<=v.MODEL['maximum_weighted_corrections']

cases=[]
with tempfile.TemporaryDirectory(prefix='phase-regression-') as tmp:
 d=Path(tmp);p=d/'haplotypes.csv'
 valid={'original_boundary':render(rows),'reverse_rows':render(rows[::-1]),'CRLF':render(rows,'\r\n'),'quoted_fields':render(rows,quoting=csv.QUOTE_ALL),'permuted_homolog_labels':render([(name,rows[-i-1][1]) for i,(name,_) in enumerate(rows)])}
 for name,data in valid.items():
  p.write_bytes(data);full(d);cases.append(dict(name=name,expected='accept',passed=True))
 bad={'leading_space':render([(rows[0][0],' '+rows[0][1])]+rows[1:]),'trailing_space':render([(rows[0][0],rows[0][1]+' ')]+rows[1:]),'unicode_digit':render([(rows[0][0],'０'+rows[0][1][1:])]+rows[1:]),'embedded_newline':render([(rows[0][0],rows[0][1][:10]+'\n'+rows[0][1][10:])]+rows[1:]),'too_short':render([(rows[0][0],rows[0][1][:-1])]+rows[1:]),'too_long':render([(rows[0][0],rows[0][1]+'0')]+rows[1:]),'missing_row':render(rows[:-1]),'duplicate_id':render(rows[:-1]+[rows[0]]),'extra_row':render(rows+[rows[0]]),'wrong_id':render([('H17',rows[0][1])]+rows[1:]),'wrong_id_space':render([(' H01',rows[0][1])]+rows[1:]),'empty':b'','extra_column':render(rows).replace(b'H01,',b'H01,x,',1),'wrong_header':render(rows).replace(b'haplotype_id',b'homolog_id',1),'blank_row':render(rows)+b'\n','BOM':b'\xef\xbb\xbf'+render(rows),'bare_CR':render(rows,'\r'),'non_binary':render([(rows[0][0],'N'+rows[0][1][1:])]+rows[1:]),'wrong_dosage':render([(rows[0][0],str(1-int(rows[0][1][0]))+rows[0][1][1:])]+rows[1:]),'too_large':b'x'*1048577,'invalid_utf8':b'\xff'}
 # Independent per-column rotations preserve every dosage but destroy read phase.
 scrambled=[''.join(truth[(h+j)%len(truth)][j] for j in range(len(truth[0]))) for h in range(len(truth))]
 v.check_dosages(scrambled);assert v.correction_cost(scrambled)>v.MODEL['maximum_weighted_corrections']
 bad['dosage_correct_bad_phase']=render([(name,scrambled[i]) for i,(name,_) in enumerate(rows)])
 for name,data in bad.items():
  p.write_bytes(data)
  try:full(d)
  except (AssertionError,OSError,ValueError,csv.Error):cases.append(dict(name=name,expected='reject',passed=True))
  else:raise AssertionError(f'Bad artifact accepted: {name}')
 p.unlink()
 for kind in ['missing','symlink','fifo','directory']:
  if kind=='symlink':
   (d/'other.csv').write_bytes(render(rows));p.symlink_to(d/'other.csv')
  elif kind=='fifo':os.mkfifo(p)
  elif kind=='directory':p.mkdir()
  try:full(d)
  except (AssertionError,OSError,ValueError,csv.Error):cases.append(dict(name=kind,expected='reject',passed=True))
  else:raise AssertionError(f'Bad artifact accepted: {kind}')
  if kind=='directory':p.rmdir()
  elif kind!='missing':p.unlink()
Path('tools/polyploid-phase/draft/regressions.json').write_text(json.dumps(dict(passed=True,cases=cases),indent=2)+'\n')
print(f'{len(valid)} valid and {len(cases)-len(valid)} invalid artifacts behaved correctly.')

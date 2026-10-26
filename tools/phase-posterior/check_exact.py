"""Independent labelled-phase enumeration checks quotient posterior inference."""
import itertools,json,math,struct,subprocess
from pathlib import Path
import numpy as np
from generate import generate,tables
DRAFT=Path('tools/phase-posterior/draft')


def enumerate_labelled(model):
 n=model['variants'];k=model['ploidy'];columns=np.zeros((20,k),np.int8)
 for i,c in enumerate(itertools.combinations(range(k),3)):columns[i,list(c)]=1
 choices=np.array(list(itertools.product(range(20),repeat=n)),np.int16)
 H=columns[choices] # labelled configuration, variant, homolog
 likelihoods=[];pzero=[];nmismatch=[]
 for r in model['reads']:
  emission=np.ones((len(H),k),np.longdouble);wrong=np.zeros((len(H),k),np.int8)
  for c in r['calls']:
   mismatch=H[:,c['variant'],:]!=c['allele'];wrong+=mismatch
   emission*=np.where(mismatch,np.longdouble(c['error_probability']),1-np.longdouble(c['error_probability']))
  l=emission.sum(1);likelihoods.append(l/k);pzero.append((emission*(wrong==0)).sum(1)/l);nmismatch.append((emission*wrong).sum(1)/l)
 mass=np.prod(likelihoods,axis=0);norm=mass.sum();mass/=norm
 answers=np.array([[np.dot(mass,p),np.dot(mass,m)] for p,m in zip(pzero,nmismatch)],float)
 return float(np.log(norm/len(H))),answers

reports=[]
for span in [2,3]:
 for seed in [71,101,208]:
  model=generate(n=4,span=span,seed=seed);tab=tables(model);stem=DRAFT/f'check-{span}-{seed}'
  with stem.with_suffix('.tables').open('wb') as f:f.write(struct.pack('<iii',tab.shape[0],2,tab.shape[-1]));f.write(tab.astype('<f8').tobytes())
  subprocess.run(['/private/tmp/phase-orbit-graph',str(span-1),str(stem.with_suffix('.graph'))],check=True,capture_output=True)
  subprocess.run(['/private/tmp/phase-infer',str(stem.with_suffix('.graph')),str(stem.with_suffix('.tables')),str(stem.with_suffix('.checkpoint')),str(stem.with_suffix('.result'))],check=True,capture_output=True)
  actual=np.fromfile(stem.with_suffix('.result'),'<f8');z,ans=enumerate_labelled(model)
  dz=abs(z-actual[0]);dp=float(np.max(abs(ans-actual[1:].reshape(-1,2))))
  assert dz<1e-11 and dp<1e-11,(span,seed,dz,dp)
  reports.append(dict(span=span,seed=seed,labelled_configurations=20**4,log_evidence_error=dz,posterior_error=dp))
  print(reports[-1],flush=True)
# Burnside independently predicts the number of unlabelled histories.
fixed=[]
for perm in itertools.permutations(range(6)):
 count=0
 for ones in itertools.combinations(range(6),3):
  s=set(ones);count+=({perm[i] for i in s}==s)
 fixed.append(count)
census={str(d):sum(f**d for f in fixed)//720 for d in range(1,9)}
assert census['7']==1822528
(DRAFT/'exhaustive-validation.json').write_text(json.dumps(dict(passed=True,cases=reports,burnside_orbit_counts=census),indent=2)+'\n')
print('Burnside counts',census)

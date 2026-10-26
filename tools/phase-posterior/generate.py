"""Original six-copy sequence experiment and model-defined emission tables."""
import csv,itertools,json,struct
from pathlib import Path
import numpy as np

def generate(n=64,span=8,seed=20261102):
 rng=np.random.default_rng(seed);H=np.zeros((6,n),np.int8)
 for j in range(n):H[rng.choice(6,3,replace=False),j]=1
 reads=[]
 for end in range(span-1,n):
  for repeat in range(2):
   h=int(rng.integers(6));js=np.arange(end-span+1,end+1);eps=rng.choice([.025,.08,.19],span)
   calls=H[h,js]^(rng.random(span)<eps)
   reads.append(dict(id=f'r{len(reads)+1:04d}',end=end,calls=[dict(variant=int(j),allele=int(a),error_probability=float(e)) for j,a,e in zip(js,calls,eps)]))
 return dict(ploidy=6,variants=n,dosages=[3]*n,span=span,reads=reads)

def tables(model):
 span=model['span'];n=model['variants'];readlist=model['reads'];tab=np.empty((n-span+1,2,3,2**span),np.float64)
 for i,r in enumerate(readlist):
  end=r['end']
  for pattern in range(2**span):
   likelihood=1.;errors=0
   for call in r['calls']:
    bit=(pattern>>(end-call['variant']))&1;wrong=bit!=call['allele'];errors+=wrong
    likelihood*=call['error_probability'] if wrong else 1-call['error_probability']
   tab[i//2,i%2,:,pattern]=[likelihood,likelihood*(errors==0),likelihood*errors]
 return tab

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--variants',type=int,default=64);p.add_argument('--span',type=int,default=8);a=p.parse_args()
 model=generate(a.variants,a.span);out=Path('tools/phase-posterior/draft')
 (out/f'model-{a.span}.json').write_text(json.dumps(model,indent=2)+'\n')
 tab=tables(model)
 with (out/f'tables-{a.span}.bin').open('wb') as f:f.write(struct.pack('<iii',tab.shape[0],2,tab.shape[-1]));f.write(tab.astype('<f8').tobytes())
 print(tab.shape)

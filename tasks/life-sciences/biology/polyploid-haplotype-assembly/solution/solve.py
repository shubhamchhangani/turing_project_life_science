"""Input-only iterated search with exact dosage-preserving moves."""
import csv,json,time
from pathlib import Path
import numpy as np
from search import costs,coordinate,anneal


def main():
    source=Path('/root/data');model=json.loads((source/'model.json').read_text())
    with (source/'variants.csv').open() as stream:variants=list(csv.DictReader(stream))
    lookup={r['variant_id']:j for j,r in enumerate(variants)}
    d=np.array([int(r['alt_dosage']) for r in variants],np.int32)
    reads={}
    with (source/'read_calls.csv').open() as stream:
        for r in csv.DictReader(stream):reads.setdefault(r['read_id'],[]).append((lookup[r['variant_id']],int(r['allele']),int(r['weight'])))
    starts=[0];positions=[];alleles=[];weights=[]
    for calls in reads.values():
        for j,a,w in calls:positions.append(j);alleles.append(a);weights.append(w)
        starts.append(len(positions))
    data=tuple(np.array(x,np.int32) for x in (starts,positions,alleles,weights))
    k=model['ploidy'];n=len(d);rng=np.random.default_rng(20261101)
    best=10**12;population=[];start=time.monotonic()
    for run in range(1500):
        if run<12 or run%10==0:
            H=np.zeros((k,n),np.int8)
            for j,v in enumerate(d):H[rng.permutation(k)[:v],j]=1
        else:
            H=population[int(rng.integers(min(6,len(population))))][1].copy()
            for move in range(int(rng.integers(2,14))):
                a,b=rng.choice(k,2,replace=False);lo=int(rng.integers(n));hi=min(n,lo+int(rng.integers(10,110)))
                H[[a,b],lo:hi]=H[[b,a],lo:hi]
            for j in rng.choice(n,size=int(rng.integers(0,30)),replace=False):rng.shuffle(H[:,j])
        H=coordinate(H,*data,d);H,score=anneal(H,*data,d,run+1,sweeps=300)
        assert int(costs(H,*data).min(1).sum())==score
        if score<best:best=score;winner=H.copy()
        population.append((score,H));population.sort(key=lambda x:x[0]);population=population[:16]
        print(json.dumps(dict(run=run,best=int(best),seconds=time.monotonic()-start)),flush=True)
        if best<=model['maximum_weighted_corrections']:break
    else:raise RuntimeError('Search did not reach the specified correction budget')
    out=Path('/root/results');out.mkdir(parents=True,exist_ok=True)
    with (out/'haplotypes.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['haplotype_id','alleles'])
        writer.writerows((f'H{i+1:02d}',''.join(map(str,h))) for i,h in enumerate(winner))

if __name__=='__main__':main()

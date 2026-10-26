"""Independent finite-field observables: vertex indexing, Householder reduction."""
import itertools
import math
import numpy as np
from scipy.linalg import inv
from pfapack.pfaffian import pfaffian
from independent import DecoratedTorus


class IndependentFieldTorus(DecoratedTorus):
    def __init__(self,width,height,sites):
        super().__init__(width,height)
        self.sites=sites
        self.primal_to_dual=[]
        for x,y,d in self.lookup:
            other=(x,(y-1)%height,1) if d==0 else ((x-1)%width,y,0)
            self.primal_to_dual.append(self.lookup[other])

    def sector_derivatives(self,weights):
        polynomial=0.
        derivative=np.zeros(len(weights))
        for tx,ty in itertools.product([0,1],repeat=2):
            matrix=self.matrix.copy()
            baseline_sign=self.match_sign
            baseline_log=0.
            entries=[]
            for weight,(i,j,sign,wx,wy) in zip(weights,self.links):
                value=sign/weight*(-1)**(tx*wx+ty*wy)
                matrix[i,j]=value;matrix[j,i]=-value
                baseline_sign*=np.sign(value)
                baseline_log+=math.log(abs(value))
                entries.append(value)
            inverse=inv(matrix,check_finite=False)
            answer=pfaffian(matrix*.5,method='H',overwrite_a=True)
            value=np.sign(answer)*baseline_sign*math.exp(math.log(abs(answer))+len(matrix)/2*math.log(2)-baseline_log)
            gradient=np.array([2*(inverse[j,i]*entry-1) for entry,(i,j,_,_,_) in zip(entries,self.links)])
            polynomial+=value/4
            derivative+=value/4*gradient
        return polynomial,derivative[self.primal_to_dual]

    def observables(self,parameters):
        k=parameters[:-3];h=parameters[-3:]
        w=np.empty(len(k))
        for index,dual in enumerate(self.primal_to_dual):w[dual]=math.exp(-2*k[index])
        masks=[0,3,5,6]
        values=[];derivatives=[]
        for mask in masks:
            selected=[p for i,p in enumerate(self.sites) if mask>>i&1][::-1]
            modified=w.copy()
            for (x,y),(xx,yy) in zip(selected[::2],selected[1::2]):
                while y!=yy:
                    y=(y-1)%self.height
                    modified[self.lookup[(x-1)%self.width,y,0]]*=-1
                while x!=xx:
                    x=(x-1)%self.width
                    modified[self.lookup[x,(y-1)%self.height,1]]*=-1
            value,gradient=self.sector_derivatives(modified)
            values.append(value);derivatives.append(gradient)
        t=np.tanh(h)
        coeff=np.array([1,t[0]*t[1],t[0]*t[2],t[1]*t[2]])
        total=coeff@values
        logz=math.log(2)+sum(k)+sum(np.log(np.cosh(h)))+math.log(total)
        edges=1+coeff@np.array(derivatives)/total
        magnetization=[]
        for i in range(3):
            change=sum(values[j]*np.prod([t[a] for a in range(3) if masks[j]>>a&1 and a!=i]) for j in range(4) if masks[j]>>i&1)
            magnetization.append(t[i]+(1-t[i]**2)*change/total)
        return logz,np.r_[edges,magnetization]


if __name__=='__main__':
    import json,time
    from pathlib import Path
    from field_inverse_prototype import FieldTorus,exhaustive
    rng=np.random.default_rng(456)
    small_checks=[]
    for w,h in [(3,3),(3,4),(4,4)]:
        sites=[(0,0),(w-1,h-1),(1,h-1)]
        solver=IndependentFieldTorus(w,h,sites)
        for lo,hi in [(.15,.35),(.35,.65),(-.2,.65)]:
            parameters=np.r_[rng.uniform(lo,hi,len(solver.links)),rng.uniform(-.6,.6,3)]
            expected=exhaustive(FieldTorus(w,h,sites),parameters)
            actual=solver.observables(parameters)
            error=max(abs(actual[0]-expected[0]),max(abs(actual[1]-expected[1])))
            small_checks.append({'width':w,'height':h,'range':[lo,hi],'max_error':float(error)})
            assert error<1e-10
    root=Path(__file__).resolve().parent/'draft/field-inverse'
    g=json.loads((root/'generation.json').read_text())
    solver=IndependentFieldTorus(g['width'],g['height'],[tuple(p) for p in g['field_sites']])
    start=time.monotonic();z,means=solver.observables(np.array(g['truth']))
    error=max(abs(means-np.array(g['observables'])))
    result={'small_checks':small_checks,'full_gradient_max_difference':float(error),'full_logz_difference':abs(z-g['log_partition']),'seconds':time.monotonic()-start,'all_passed':error<1e-10 and abs(z-g['log_partition'])<1e-10}
    result['all_passed']=bool(result['all_passed'])
    (root/'independent-validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result,flush=True)

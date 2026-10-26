"""Cross-check algorithms and analytic photon/zero-transcription limits."""
import copy,json,math
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import poisson
import generate,reference,independent

def run():
    rng=np.random.default_rng(20261026);cases=[]
    for k in range(8):
        m=generate.model(2+k%2,2+(k//2)%2);m['initial']=[k%2,1,1];m['count_limit']=5
        for branch in ['mother','daughter_a','daughter_b']:
            m[branch]=[dict(duration=float(rng.uniform(.15,.65)),stimulus=float(rng.uniform(.3,1.8)),photons=int(rng.integers(0,4))) for _ in range(2)]
        a=reference.solve(m);b=independent.solve(m)
        err=max(float(np.max(np.abs(np.array(a[key])-np.array(b[key])))) for key in ['parent','daughter_rna','counts','log_evidence','count_overflow'])
        assert err<2e-11
        assert abs(sum(a['parent'])-1)<1e-12 and abs(np.sum(a['daughter_rna'])-1)<1e-12
        cases.append({'case':k,'max_difference':err})
    m=generate.model(3,4);m['initial']=[1,1,2];m['rates']['photon_brightness']=0;m['rates']['transcription']=0
    m['count_limit']=3;a=reference.solve(m);b=independent.solve(m)
    expected=sum(poisson.logpmf(t['photons'],m['rates']['photon_background']*t['duration']) for name in ['mother','daughter_a','daughter_b'] for t in m[name])
    for result in [a,b]:
        assert abs(result['log_evidence']-expected)<2e-11
        assert abs(result['counts'][0]-1)<2e-12 and max(abs(v) for v in result['counts'][1:])==0
        assert abs(result['count_overflow'])<2e-12
    # State-independent photon observations must leave the mother prior unchanged.
    st=reference.states(m);prior=np.zeros(len(st));prior[st.index(tuple(m['initial']))]=1
    for t in m['mother']:
        Q,_,_=reference.generator(m,t['stimulus']);prior=expm_multiply(Q*t['duration'],prior)
    assert np.max(np.abs(prior-a['parent']))<2e-12
    out={'passed':True,'cross_algorithm_cases':cases,'analytic_limits_passed':True}
    (Path(__file__).parent/'draft').mkdir(exist_ok=True)
    (Path(__file__).parent/'draft/small-validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)

if __name__=='__main__':run()

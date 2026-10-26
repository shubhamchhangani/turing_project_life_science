"""Validate image integration and independent photon-time likelihood calculations."""
import json,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.spatial.transform import Rotation
from model import projection,photon_model
from likelihood import load,prepare,trace_loglike

def direct(Q,pi,rates,times,channels,duration):
    state=pi.copy();logp=0.;previous=0.;A=Q-np.diag(rates.sum(axis=0))
    for t,c in zip(times,channels):
        state=expm(A*(t-previous))@state;state*=rates[c]
        scale=state.sum();state/=scale;logp+=np.log(scale);previous=t
    return float(logp+np.log(np.sum(expm(A*(duration-previous))@state)))

def uniformized(Q,pi,rates,times,channels,duration):
    A=Q-np.diag(rates.sum(axis=0));rate=max(-A.diagonal());P=np.eye(len(pi))+A/rate
    def advance(v,dt):
        pieces=max(1,int(np.ceil(rate*dt/12)));mu=rate*dt/pieces
        for _ in range(pieces):
            term=v.copy();w=np.exp(-mu);out=w*term
            for k in range(1,160):
                term=P@term;w*=mu/k;out+=w*term
                if k>mu+20 and w<1e-20:break
            v=out
        return v
    v=pi.copy();ll=0.;previous=0.
    for t,c in zip(times,channels):
        v=advance(v,t-previous)*rates[c];z=v.sum();ll+=np.log(z);v/=z;previous=t
    return float(ll+np.log(advance(v,duration-previous).sum()))

if __name__=='__main__':
    started=time.monotonic();root=Path(__file__).parent/'draft';c,traces=load(root/'public');rng=np.random.default_rng(831)
    # These arbitrary validation parameters are independent of the data generator.
    angles=rng.uniform([-.9,-.6],[1.7,.8],(4,2));energies=np.array([.1,.3,.7]);barriers=np.arange(6)*.2+3.1
    cases=[]
    for meta,ts,cs in traces[:8]:
        ts=ts[:12];cs=cs[:12];duration=float(ts[-1]+2e-5)
        Q,pi,rates=photon_model(c,angles,energies,barriers,meta['temperature'],meta['label_pair'])
        values=[trace_loglike(ts,cs,duration,*prepare(Q,pi,rates)),direct(Q,pi,rates,ts,cs,duration),uniformized(Q,pi,rates,ts,cs,duration)]
        error=float(np.ptp(values));assert error<2e-10;cases.append(error)
    # Zero photons: the final no-photon survival interval must still be included.
    Q,pi,rates=photon_model(c,angles,energies,barriers,300.,0)
    empty=trace_loglike(np.array([]),np.array([],dtype=np.int8),1e-4,*prepare(Q,pi,rates))
    assert abs(empty-direct(Q,pi,rates,[],[],1e-4))<1e-10 and empty<0
    # Independent direct z integration of the 3D normalized Gaussian mixture.
    test=dict(c,transfer_decay=0.);xyz=rng.normal(size=(9,3))*8;rot=Rotation.random(random_state=rng).as_matrix();shift=[1.2,-.7]
    image=projection(test,xyz,rot,shift,0.);xyz=xyz@rot.T;axis=(np.arange(c['image_size'])-(c['image_size']-1)/2)*c['pixel_size']
    z=np.linspace(-100,100,20001);sigma=c['gaussian_sigma'];errors=[]
    for iy,ix in rng.integers(7,25,(12,2)):
        density=np.zeros(len(z))
        for pos,weight in zip(xyz,c['bead_weights']):
            d2=(axis[ix]-pos[0]-shift[0])**2+(axis[iy]-pos[1]-shift[1])**2+(z-pos[2])**2
            density+=weight*np.exp(-d2/(2*sigma**2))/(2*np.pi*sigma**2)**1.5
        errors.append(abs(np.trapz(density,z)-image[iy,ix]))
    assert max(errors)<1e-12
    result=dict(passed=True,photon_likelihood_max_log_difference=max(cases),zero_photon_tail_checked=True,
        projection_quadrature_max_difference=float(max(errors)),seconds=time.monotonic()-started,
        scope='forward numerical checks only; no inverse validation or difficulty claim')
    (root/'forward-validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

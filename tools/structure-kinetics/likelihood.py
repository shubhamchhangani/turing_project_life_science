"""Photon-time likelihood and public-data loader. No private truth access."""
import json
from pathlib import Path
import numpy as np
from numba import njit
from scipy.linalg import eigh
from model import photon_model

@njit(cache=True)
def trace_loglike(times,channels,duration,eigenvalues,transforms,start):
    a=start.copy();previous=0.;ll=0.
    for i in range(len(times)):
        a*=np.exp(eigenvalues*(times[i]-previous))
        a=transforms[channels[i]]@a
        normalizer=np.dot(start,a)
        if normalizer<=0:return -np.inf
        a/=normalizer;ll+=np.log(normalizer);previous=times[i]
    tail=np.dot(start,np.exp(eigenvalues*(duration-previous))*a)
    if tail<=0:return -np.inf
    return ll+np.log(tail)

def prepare(Q,pi,rates):
    root=np.sqrt(pi)
    # alpha = diag(sqrt(pi)) beta gives a symmetric reversible killed generator.
    symmetric=Q*root[None,:]/root[:,None]-np.diag(rates.sum(axis=0))
    assert np.max(np.abs(symmetric-symmetric.T))<1e-6
    values,vectors=eigh(symmetric)
    transforms=np.stack([vectors.T@(r[:,None]*vectors) for r in rates])
    return values,transforms,np.ascontiguousarray(vectors.T@root)

def load(public):
    public=Path(public);c=json.loads((public/'config.json').read_text());meta=json.loads((public/'traces.json').read_text())
    with np.load(public/'photons.npz',allow_pickle=False) as photons:
        traces=[]
        for row in meta:
            mask=photons['trace_id']==row['trace_id']
            traces.append((row,photons['time'][mask].copy(),photons['channel'][mask].copy()))
    return c,traces

def decode(theta):
    return np.asarray(theta[:8]).reshape(4,2),np.asarray(theta[8:11]),np.asarray(theta[11:17])

def loglike(config,traces,theta):
    angles,energies,barriers=decode(theta);cache={};value=0.
    for meta,times,channels in traces:
        key=(meta['temperature'],meta['label_pair'])
        if key not in cache:cache[key]=prepare(*photon_model(config,angles,energies,barriers,*key))
        value+=trace_loglike(times,channels,meta['duration'],*cache[key])
    return float(value)

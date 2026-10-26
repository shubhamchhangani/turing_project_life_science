"""Independent compiled sequential-update implementation for provenance."""
import numpy as np
from numba import njit


@njit
def source_moments(inverse, paths, delta, weights, count):
    out=np.zeros(count);out[0]=1;size=len(weights)
    for mask in range(1,count):
        parity=0;word0=np.uint64(0);word1=np.uint64(0);bits=mask;bit=0
        while bits:
            if bits&1:parity^=1;word0^=paths[bit,0];word1^=paths[bit,1]
            bits>>=1;bit+=1
        if parity:continue
        chosen=np.empty(size,np.int64);k=0
        for edge in range(size):
            present=(word0>>np.uint64(edge))&np.uint64(1) if edge<64 else (word1>>np.uint64(edge-64))&np.uint64(1)
            if present:chosen[k]=edge;k+=1
        if not k:out[mask]=1;continue
        n=2*k;b=np.empty((n,n),np.complex128)
        for i in range(n):
            for j in range(n):b[i,j]=inverse[2*chosen[i//2]+i%2,2*chosen[j//2]+j%2]
        ratio=1+0j
        for pair in range(k):
            p=2*pair;edge=chosen[pair];d=delta[edge];q=1-d*b[p,p+1]
            ratio*=q/weights[edge];scale=d/q
            for i in range(p+2,n):
                for j in range(i+1,n):
                    b[i,j]+=scale*(b[i,p]*b[j,p+1]-b[i,p+1]*b[j,p]);b[j,i]=-b[i,j]
        out[mask]=ratio.real
    return out


@njit
def source_invert(values):
    n=len(values);stride=n//2
    while stride:
        for start in range(0,n,2*stride):
            for j in range(stride):
                a=values[start+j];b=values[start+j+stride]
                values[start+j]=a+b;values[start+j+stride]=a-b
        stride//=2
    values/=n

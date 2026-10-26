"""Compiled, bounded-memory signed moment enumeration and exact reweighting."""
import numpy as np
from numba import njit


@njit
def moments(inverse, paths, delta, factor, count):
    out=np.zeros(count);out[0]=1
    width=len(delta)
    for mask in range(1,count):
        bits=mask;parity=0;lo=np.uint64(0);hi=np.uint64(0);bit=0
        while bits:
            if bits&1:
                parity^=1;lo^=paths[bit,0];hi^=paths[bit,1]
            bits>>=1;bit+=1
        if parity:continue
        chosen=np.empty(width,np.int64);k=0
        for e in range(width):
            yes=(lo>>np.uint64(e))&np.uint64(1) if e<64 else (hi>>np.uint64(e-64))&np.uint64(1)
            if yes:chosen[k]=e;k+=1
        if k==0:out[mask]=1;continue
        n=2*k;b=np.empty((n,n),np.complex128);value=complex(-1 if k%2 else 1)
        for i in range(n):
            for j in range(n):b[i,j]=inverse[2*chosen[i//2]+i%2,2*chosen[j//2]+j%2]
        for i in range(k):
            e=chosen[i];b[2*i,2*i+1]-=1/delta[e];b[2*i+1,2*i]+=1/delta[e];value*=factor[e]
        for p in range(0,n-1,2):
            pivot=p+1
            for j in range(p+2,n):
                if abs(b[p,j])>abs(b[p,pivot]):pivot=j
            if b[p,pivot]==0:value=0j;break
            if pivot!=p+1:
                for j in range(n):b[p+1,j],b[pivot,j]=b[pivot,j],b[p+1,j]
                for i in range(n):b[i,p+1],b[i,pivot]=b[i,pivot],b[i,p+1]
                value=-value
            d=b[p,p+1];value*=d
            for i in range(p+2,n):
                for j in range(i+1,n):
                    b[i,j]+=(b[p+1,i]*b[p,j]-b[p,i]*b[p+1,j])/d
                    b[j,i]=-b[i,j]
        out[mask]=value.real
    return out


@njit
def character_inverse(values):
    step=1;n=len(values)
    while step<n:
        for start in range(0,n,2*step):
            for j in range(start,start+step):
                a=values[j];b=values[j+step];values[j]=a+b;values[j+step]=a-b
        step*=2
    values/=n


@njit
def reweight(values, field_sites, field_values, edge_a, edge_b, coupling, nlig, clamp, nreport):
    out=np.zeros(1<<nreport);total=0.;clamp_total=0.
    for mask in range(len(values)):
        energy=0.
        for j in range(len(field_sites)):energy+=field_values[j]*(1-2*((mask>>field_sites[j])&1))
        for j in range(len(coupling)):
            energy+=coupling[j]*(1-2*((mask>>edge_a[j])&1))*(1-2*((mask>>edge_b[j])&1))
        value=values[mask]*np.exp(energy);total+=value
        if mask&((1<<nlig)-1)==clamp:
            clamp_total+=value;out[(mask>>nlig)&((1<<nreport)-1)]+=value
    return total,clamp_total,out/clamp_total


def enumerate_moments(engine):
    assert len(engine.union)<=128
    paths=np.zeros((len(engine.sites),2),np.uint64)
    for i,path in enumerate(engine.paths):
        for edge in path:
            e=engine.lookup[edge];paths[i,e//64]^=np.uint64(1)<<np.uint64(e%64)
    edges=np.array(engine.union);delta=engine.r[2*edges]*(engine.t[edges]-1/engine.t[edges])
    return moments(engine.inverse,paths,delta,delta/engine.t[edges],1<<len(engine.sites))

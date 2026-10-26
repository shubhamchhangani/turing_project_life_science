"""Signed spin moments using skew Kac-Ward matrices and Pfaffian updates."""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from engine import Planar


def pfaffian(matrix):
    a=matrix.copy();n=len(a);result=1+0j
    for k in range(0,n-1,2):
        pivot=k+1+int(np.argmax(np.abs(a[k,k+1:])))
        if a[k,pivot]==0:return 0j
        if pivot!=k+1:
            a[[k+1,pivot],:]=a[[pivot,k+1],:]
            a[:,[k+1,pivot]]=a[:,[pivot,k+1]]
            result=-result
        value=a[k,k+1];result*=value
        if k+2<n:
            u=a[k,k+2:].copy();v=a[k+1,k+2:].copy()
            a[k+2:,k+2:]+=(np.outer(v,u)-np.outer(u,v))/value
    return result


class Signed(Planar):
    def __init__(self,model):
        super().__init__(model)
        xy=np.array(model['coordinates'])
        vectors=(xy[self.dst]-xy[self.src])@np.array([1,1j])
        self.r=vectors/np.abs(vectors)
        d=np.arange(self.dim)
        a=csc_matrix((np.r_[self.phases*self.r[self.rows^1],self.r[d^1]/self.t[d//2]],
                      (np.r_[self.rows^1,d^1],np.r_[self.cols,d])),shape=(self.dim,self.dim))
        assert abs(a+a.T).max()<1e-12
        self.sites=model['ligand_sites']+model['reporter_sites']
        self.paths=[set(self.path(self.sites[0],site)) for site in self.sites]
        self.union=sorted(set().union(*self.paths))
        selected=np.array([d for e in self.union for d in [2*e,2*e+1]])
        rhs=np.zeros((self.dim,len(selected)),complex);rhs[selected,np.arange(len(selected))]=1
        self.inverse=splu(a).solve(rhs)[selected,:]
        self.inverse=(self.inverse-self.inverse.T)/2
        self.lookup={e:i for i,e in enumerate(self.union)}

    def moment(self,mask):
        if mask.bit_count()%2:return 0.0
        changed=set()
        for j,path in enumerate(self.paths):
            if mask>>j&1:changed.symmetric_difference_update(path)
        if not changed:return 1.0
        edges=sorted(changed);k=len(edges)
        local=np.array([d for e in edges for d in [2*self.lookup[e],2*self.lookup[e]+1]])
        b=self.inverse[np.ix_(local,local)].copy()
        delta=self.r[2*np.array(edges)]*(self.t[edges]-1/self.t[edges])
        j=np.arange(k)*2
        b[j,j+1]-=1/delta;b[j+1,j]+=1/delta
        value=(-1)**k*pfaffian(b)*np.prod(delta/self.t[edges])
        assert abs(value.imag)<1e-8
        assert abs(value.real)<1+1e-8
        return float(value.real)


def distribution(model,progress=None):
    engine=Signed(model);count=1<<len(engine.sites)
    values=np.empty(count)
    for mask in range(count):
        values[mask]=engine.moment(mask)
        if progress and mask%1024==0:progress(mask,count)
    step=1
    while step<count:
        blocks=values.reshape(-1,2*step);left=blocks[:,:step].copy();right=blocks[:,step:].copy()
        blocks[:,:step]=left+right;blocks[:,step:]=left-right;step*=2
    values/=count
    nlig=len(model['ligand_sites']);clamp=sum((s==-1)<<i for i,s in enumerate(model['ligand_spins']))
    joint=values[clamp::1<<nlig];pc=joint.sum()
    return {'log_partition':engine.logz,'log_clamp_probability':float(np.log(pc)),
            'probabilities':(joint/pc).tolist(),'minimum_joint_probability':float(values.min())}

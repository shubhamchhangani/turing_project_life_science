import numpy as np, pandas as pd, sys, pysam
from scipy.stats import binom
from scipy.optimize import minimize
p=float(sys.argv[1]) if len(sys.argv)>1 else 0.607; f=float(sys.argv[2]) if len(sys.argv)>2 else 0.205
et=1.0e-4; en=1.65e-4
df=pd.read_pickle('/root/work/post.pkl')
fa=pysam.FastaFile('/root/data/reference.fa')
seqs={c:fa.fetch(c).upper() for c in fa.references}
def cpg(c,pos,r):
    if c=='pseudo': return False
    s=seqs[c]
    if r=='C': return pos<len(s) and s[pos]=='G'
    if r=='G': return pos>=2 and s[pos-2]=='C'
    return False
df['cpg']=[cpg(c,pp,r) for c,pp,r in zip(df.chrom,df.pos,df.ref)]
# pseudo rows: split by approximate CpG fraction of C/G sites
ps=df.chrom=='pseudo'
cg=df[~ps&df.ref.isin(['C','G'])]
# CpG fraction among C/G reference positions
tot=sum(s.count('C')+s.count('G') for s in seqs.values()); ncpg=sum(2*s.count('CG') for s in seqs.values()); fr=ncpg/tot
print('cpg frac',fr)
extra=df[ps&df.ref.isin(['C','G'])].copy(); extra['cpg']=True; extra['w']=extra.w*fr
df.loc[ps&df.ref.isin(['C','G']),'w']*= (1-fr)
df=pd.concat([df,extra],ignore_index=True)
df['k']=df.cl+np.where(df.cpg,'_CpG','')
W=df.w.values; ac=df.ac.values.astype(float); Dc=df.Dc.values.astype(float); an=df.an.values.astype(float); Dn=df.Dn.values.astype(float); C=df.C.values.astype(float)
kv=df.k.values; ks=sorted(set(kv)); kidx=np.array([ks.index(x) for x in kv])
comps=[('mal',0.25,0.25),('mal',0.55,0.55),('mal',1,1),('mal',1.5,1.5),('mal',2,2),('mal',3,3),('non',0.3,1.1),('non',0.6,1.8),('non',0.13,0.55)]
ismal=np.array([c[0]=='mal' for c in comps])
def lik(mt,mn,mal):
    vt=p*mt/(2*(1-p)+p*C); vn=f*mn/(2*(1-f)+f*C)
    l=binom.logpmf(ac,Dc,np.minimum(vt+et,1))+binom.logpmf(an,Dn,np.minimum(vn+en,1))
    if mal: l=np.where(mt<=C+1e-9,l,-np.inf)
    return l
L0=binom.logpmf(ac,Dc,et)+binom.logpmf(an,Dn,en)
w=np.ones(len(comps)); w[ismal]/=ismal.sum(); w[~ismal]/=(~ismal).sum()
pm=np.full(len(ks),2e-4); pn=np.full(len(ks),1e-4)
for it in range(80):
    L=np.stack([lik(a,b,t=='mal') for t,a,b in comps],1)
    pri=np.where(ismal[None,:],pm[kidx][:,None],pn[kidx][:,None])*w[None,:]
    A=np.log(pri)+L; B=np.log1p(-pm[kidx]-pn[kidx])+L0
    mx=np.maximum(A.max(1),B); R=np.exp(A-mx[:,None]); r0=np.exp(B-mx); R/=(R.sum(1)+r0)[:,None]
    Rm=R[:,ismal].sum(1); Rn=R[:,~ismal].sum(1)
    for j in range(len(ks)):
        m_=kidx==j; pm[j]=(Rm[m_]*W[m_]).sum()/W[m_].sum(); pn[j]=(Rn[m_]*W[m_]).sum()/W[m_].sum()
    ww=(R*W[:,None]).sum(0); w=ww.copy(); w[ismal]/=ww[ismal].sum(); w[~ismal]/=ww[~ismal].sum()
    new=[]
    for j,(t,a,b) in enumerate(comps):
        r=R[:,j]*W; sel=r>1e-4
        if t=='mal' and a in (1,2,3):  new.append((t,a,b)); continue
        acs,Dcs,ans,Dns,Cs,rs=ac[sel],Dc[sel],an[sel],Dn[sel],C[sel],r[sel]
        def nl(x):
            mt=x[0]; mn=x[0] if t=='mal' else x[1]
            vt=p*mt/(2*(1-p)+p*Cs); vn=f*mn/(2*(1-f)+f*Cs)
            return -(rs*(binom.logpmf(acs,Dcs,np.clip(vt+et,1e-9,1))+binom.logpmf(ans,Dns,np.clip(vn+en,1e-9,1)))).sum()
        x0=[a] if t=='mal' else [a,b]
        x=minimize(nl,x0,method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-3)).x
        new.append((t,x[0],x[0] if t=='mal' else x[1]))
    comps=new
print('p f',p,f)
for (t,a,b),wv,tw in zip(comps,w,ww): print(t,round(a,3),round(b,3),round(tw,1))
for j,k in enumerate(ks): print(k,'mal n',round((R[kidx==j][:,ismal].sum(1)*W[kidx==j]).sum(),1),'non n',round((R[kidx==j][:,~ismal].sum(1)*W[kidx==j]).sum(),1))
df['pmal']=R[:,ismal].sum(1); df['pnon']=R[:,~ismal].sum(1)
ll=(np.log(np.exp(A-mx[:,None]).sum(1)+np.exp(B-mx))+mx); print('LL',(ll*W).sum())
totm=(df.pmal*W).sum()
d=df[df.chrom!='pseudo']; d.to_pickle(f'/root/work/post4_{p}_{f}.pkl')
print('expected mal',totm,'non',(df.pnon*W).sum())
for th in [0.3,0.4,0.5,0.6,0.7,0.9]:
    s=d[d.pmal>th]; print(th,len(s),'expTP',round(s.pmal.sum(),1),'prec',round(s.pmal.mean(),3),'recall',round(s.pmal.sum()/totm,3))

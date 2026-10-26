import pickle,numpy as np,pandas as pd,pysam
N=pickle.load(open('/root/work/normal.scan2.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan2.pkl','rb'))
ref=pysam.FastaFile('/root/data/reference.fa')
rows=[]
z=np.zeros(6,int)
for c in T:
    seq=ref.fetch(c).upper()
    dt,at=T[c]; dn,an=N[c]
    keys=set(k for k,v in at.items() if v[0]+v[1]+v[3]>=2)|set(k for k,v in an.items() if v[0]+v[1]+v[3]>=2)
    for (p,a) in keys:
        v=at.get((p,a),z); nv=an.get((p,a),z)
        rows.append((c,p+1,seq[p],a,*v[:5],dt[p],*nv[:5],dn[p]))
df=pd.DataFrame(rows,columns=['chrom','pos','ref','alt','tF1','tF2','tlow','tdisc','tconc','dt','nF1','nF2','nlow','ndisc','nconc','dn'])
df['talt']=df.tF1+df.tF2; df['nalt']=df.nF1+df.nF2
df.to_csv('/root/work/cand2.tsv',sep='\t',index=False)
print(len(df))

import pickle,numpy as np,pandas as pd,pysam
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan.pkl','rb'))
ref=pysam.FastaFile('/root/data/reference.fa')
rows=[]
for c in T:
    seq=ref.fetch(c).upper()
    dt,at=T[c]; dn,an=N[c]
    for (p,a),v in at.items():
        t=v[0]+v[1]
        n=an.get((p,a)); nv=n if n is not None else np.zeros(4,int)
        if t>=2 or (nv[0]+nv[1])>=2:
            rows.append((c,p+1,seq[p],a,v[0],v[1],v[2],dt[p],nv[0]+nv[1],nv[2],dn[p]))
df=pd.DataFrame(rows,columns=['chrom','pos','ref','alt','tF1','tF2','tlow','dt','nalt','nlow','dn'])
df['talt']=df.tF1+df.tF2
df.to_csv('/root/work/cand.tsv',sep='\t',index=False)
print(len(df))

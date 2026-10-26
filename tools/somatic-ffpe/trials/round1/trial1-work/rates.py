import pickle,numpy as np,pysam,collections
ref=pysam.FastaFile('/root/data/reference.fa')
for s in ['tumor','normal']:
    S=pickle.load(open(f'/root/work/{s}.scan2.pkl','rb'))
    obs=collections.Counter(); dep=collections.Counter()
    hist=collections.Counter()
    for c in S:
        seq=np.frombuffer(ref.fetch(c).upper().encode(),dtype='S1')
        d,a=S[c]
        for b in 'ACGT': dep[b]+=int(d[seq==b.encode()].sum())
        for (p,alt),v in a.items():
            if v[0]+v[1] > 0.1*d[p]: continue  # skip real variants
            r=seq[p].decode()
            obs[(r,alt,'F1')]+=int(v[0]); obs[(r,alt,'F2')]+=int(v[1])
            if r=='C' and alt=='T': hist[int(v[0])]+=1
    print(s)
    for k in sorted(obs): print('  ',k,obs[k],'%.2e'%(obs[k]/(dep[k[0]]/2)))
    print(' C>T F1 per-site hist',sorted(hist.items())[:15])

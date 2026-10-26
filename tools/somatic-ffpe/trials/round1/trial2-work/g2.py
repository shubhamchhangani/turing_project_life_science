import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
c=pd.read_pickle('/root/work/cf4.pkl'); x=c[c.call]
print(np.histogram(x.dll.clip(-60,60),bins=24))
print(x[(x.dll.abs()<8)][['chrom','pos','seg','sub','tumor_afr','tumor_nfr','tv','normal_afr','normal_nfr','nv','anti','dll']].sort_values('dll').round(3).to_string())

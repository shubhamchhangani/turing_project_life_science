import pandas as pd, numpy as np
pd.set_option('display.width',250)
df=pd.read_pickle('/root/work/all.pkl')
df=df[df.alt!='.']
x=df[(df.talt<3)&(df.nalt>=2)&(df.talt>=1)]
print(x.to_string())

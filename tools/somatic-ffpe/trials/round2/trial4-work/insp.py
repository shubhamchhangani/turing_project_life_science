import numpy as np,sys
d=np.load('/root/work/calls.npz'); lo,hi=float(sys.argv[1]),float(sys.argv[2])
m=(d['llr']>=lo)&(d['llr']<hi)&(~d['germ'])
for i in np.where(m)[0]:
    print(d['cand'][i],round(d['llr'][i],1),'sa/sd',d['sa'][i],d['sd'][i],'an/N',d['an'][i],d['N'][i],'ct',int(d['isct'][i]),'cn',int(d['cn'][i]),'v',round(d['bv'][i],3))

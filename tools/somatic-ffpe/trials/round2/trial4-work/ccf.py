import numpy as np
d=np.load('/root/work/calls.npz')
sa,sd,an,N,ct,cn,germ,llr=d['sa'],d['sd'],d['an'],d['N'],d['isct'],d['cn'],d['germ'],d['llr']
p,f=0.615,0.2
DT=p*cn+2*(1-p); DN=f*cn+2*(1-f)
# MLE of x = m*ccf: tumour expects x*p/DT, normal x*f/DN ; approximate via grid
xs=np.linspace(0.01,3,300)
m=(~germ)&(llr>12)
L=np.array([ (sa[m]*np.log(x*p/DT[m])+(sd[m]-sa[m])*np.log(1-x*p/DT[m]+1e-9)+an[m]*np.log(x*f/DN[m])+(N[m]-an[m])*np.log(1-x*f/DN[m]+1e-9)) for x in xs])
xh=xs[L.argmax(0)]
print(m.sum())
h=np.histogram(xh,bins=np.arange(0,2.6,0.05)); 
for a,b in zip(h[1],h[0]): print(round(a,2),b,'#'*b)
cnm=cn[m]
for c in [2,3,4,5]:
    print(c, np.histogram(xh[cnm==c],bins=np.arange(0,3.2,0.1))[0])

import sys,os
import numpy as np

if len(sys.argv)<2:
    print 'Usage: datfile [norm]'
    print 'The last col is y'
    exit(1)

eps=1e-8

norm=0
if len(sys.argv)>2:
    norm=int(sys.argv[2])

x=[]
y=[]
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    x.append([float(m) for m in temp[0:-1]])
    y.append(float(temp[-1]))
fi.close()

xarr=np.array(x)

if norm:
    ncol=xarr.shape[1]
    for i in range(ncol):
	factor=max(xarr[:,i])-min(xarr[:,i])
	if factor<eps:
	    factor=xarr[0,i]
	xarr[:,i]/=factor

coef,residuals,rank,s = np.linalg.lstsq(xarr, y)

for i,x in enumerate(coef):
    print i,x

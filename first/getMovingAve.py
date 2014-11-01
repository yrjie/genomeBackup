#!/home/ruijie/bin/python

import sys,os

if len(sys.argv)<2:
    print 'Usage: pythont getMovingAve.py infile'
    exit(1)

def movingAve(x,n):
    ret=[]
    n/=2
    for i in xrange(len(x)+1):
    	a=max(i-n,0)
    	b=min(i+n,len(x)-1)
    	ret.append(1.0*sum(x[a:b+1])/(b-a+1))
    return ret

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(i) for i in line.split()]
    res=movingAve(temp,5)
    print "\t".join([str(i) for i in res])
fi.close()

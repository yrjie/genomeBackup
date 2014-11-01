import scipy.stats as ss
import sys,os

if len(sys.argv)<3:
    print 'Usage: datfile distr'
    exit(0)

distr=getattr(ss,sys.argv[2])
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    dat=[float(x) for x in temp]
    par=distr.fit(dat)
    print [distr.sf(x,par[0],par[1]) for x in dat]
fi.close()

import sys,os
import numpy as np

if len(sys.argv)<3:
    print 'Usage: datfile perc(0-100)'
    exit(1)

dat=[]
fi=open(sys.argv[1])
for line in fi:
    line =line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    dat.append(float(temp[0]))
fi.close()

print len(dat),np.percentile(dat, float(sys.argv[2]))

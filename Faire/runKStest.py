import sys,os
from scipy.stats import ks_2samp
import numpy as np

if len(sys.argv)<2:
    print 'Usage: datfile'
    exit(1)

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(x) for x in line.split('\t')]
#    test=[np.mean(temp)]*len(temp)
    test=[np.mean(temp)]*3
    print ks_2samp(temp, test)[1]
#    print stats.kstest(temp, 'uniform')
fi.close()


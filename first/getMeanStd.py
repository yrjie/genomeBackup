#!/home/ruijie/bin/python

import sys, os
import numpy as np

if len(sys.argv)<2:
    print 'Usage: python getMeanStd.py infile'
    exit(1)

all=[float(x.strip()) for x in open(sys.argv[1],'r').readlines()]
a=np.mean(all)
b=np.std(all)
all=[x for x in all if abs((x-a)/b)<2]
a=np.mean(all)
b=np.std(all)
print str(len(all))+"\t"+str(a)+"\t"+str(b)

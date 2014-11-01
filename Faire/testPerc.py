import sys,os
import numpy as np

if len(sys.argv)<2:
    print 'Usage: num'
    exit(1)

num=int(sys.argv[1])
a=[-i for i in range(num)]
b=0
a.sort()
for i in range(100):
    b+=a[i*num/100-1]

#for i in range(1):
#    b+=np.percentile(a, i)

import sys,os
import resource

if len(sys.argv)<2:
    print 'Usage: numberOfEntries'
    exit(1)

num=int(sys.argv[1])
a={}
one=num/10
for i in range(num):
    if i%one==0:
    	print '\t'.join([str(i),str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)])
    a[str(i)*30]=[i, False]

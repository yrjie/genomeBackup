import sys,os
import numpy

if len(sys.argv)<2:
    print 'Usage: tabFile'
    exit(1)

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    temp=[float(x) for x in line.split()]
    print numpy.std(temp)/(numpy.mean(temp)+0.1)
fi.close()

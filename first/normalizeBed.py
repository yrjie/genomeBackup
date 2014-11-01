#!/usr/bin/python

import os,sys

if len(sys.argv)<3:
    print 'Usage: normalizeBed bedfile colOfSig scale'
    exit(1)

fi=open(sys.argv[1])
lines=[]
col=int(sys.argv[2])-1
scale=int(sys.argv[3])
sum=0
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    lines.append(temp)
    sum+=float(temp[col])

for x in lines:
    x[col]=str(float(x[col])/sum*scale)
    print '\t'.join(x)


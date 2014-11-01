import scipy.stats as ss
import sys,os

if len(sys.argv)<3:
    print 'Usage: overlapCnt chromHmmCnt'
    exit(1)

overlap={}
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    overlap[temp[1]]=int(temp[0])
fi.close()

fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    if temp[1] not in overlap:
    	rate=0
    else:
    	rate=1.0*overlap[temp[1]]/int(temp[0])
    print '\t'.join([temp[1],str(rate)])
fi.close()

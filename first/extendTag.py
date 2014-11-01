#!/home/ruijie/bin/python

import sys,os
import random
import string

if len(sys.argv)<4:
    print 'Usage: python extendTag.py bamFile genomeFile extLen'
    exit(1)

mapChr={}
fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    mapChr[temp[0]]=int(temp[1])
fi.close()
extLen=int(sys.argv[3])
bedF='/tmp/'+''.join(random.choice(string.letters) for i in xrange(6))+'.bed'
cmd='bamToBed -i '+sys.argv[1]+' >'+bedF
os.system(cmd)
fi=open(bedF)
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[5]=='+':
    	temp[2]=str(min(int(temp[2])+extLen,mapChr[temp[0]]))
    else:
    	temp[1]=str(max(int(temp[1])-extLen,1))
    print '\t'.join(temp)
fi.close()
os.system('rm '+bedF)

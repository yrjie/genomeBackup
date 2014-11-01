import sys,os
import re

if len(sys.argv)<2:
    print 'Usage: smpAnnotFile'
    print '  1st col is the id, 2 row is the html'
    exit(1)

annot=[]
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    m=re.search('group: (.*?)<br>?', temp[1])
    if m==None:
    	continue
    annot.append(m.group(1)+'_'+temp[0])
fi.close()
print '\t'.join(annot)

import sys,os
import re

if len(sys.argv)<3:
    print 'Usage: idValFile annot [geneSymFile]'
    print '  annot: id symb(terms=xxx)'
    exit(1)

outfile=''
if len(sys.argv)>3:
    outfile=sys.argv[3]

mapG={}

fi=open(sys.argv[2])
for line in fi:
    #line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    m=re.search('terms=(.*?\")', temp[1].strip())
    sym=m.group(1).replace('"','')
    mapG[temp[0]]=sym
fi.close

fi=open(sys.argv[1])
fo=None
cnt=1
if outfile!='':
    fo=open(outfile,'w')
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[0] not in mapG:
    	continue
#    print '\t'.join([mapG[temp[0]],temp[1]])
    if fo:
    	#fo.write(mapG[temp[0]]+"."+str(cnt)+'\n')
    	fo.write(mapG[temp[0]]+'\n')
    	cnt+=1
    else:
    	print temp[1]
fi.close()
if fo:
    fo.close()

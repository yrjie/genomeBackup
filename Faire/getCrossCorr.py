import sys,os
import numpy as np

if len(sys.argv)<2:
    print 'Usage: aln.bed'
    exit(1)

fi=open(sys.argv[1])
mapP={}
mapM={}
chr=''
cnt={}
total=0
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if chr!=temp[0]:
    	chr=temp[0]
    	mapP[chr]={}
    	mapM[chr]={}
    	cnt[chr]=0
    if temp[5]=='+':
    	beg=int(temp[1])
	if beg not in mapP[chr]:
	    mapP[chr][beg]=1
	else:
	    mapP[chr][beg]+=1
    else:
    	beg=int(temp[2])
	if beg not in mapM[chr]:
	    mapM[chr][beg]=1
	else:
	    mapM[chr][beg]+=1
    cnt[chr]+=1
    total+=1
fi.close()

for flen in range(-50,301,10):
    a=[]
    b=[]
    ans=0
    for chr in cnt:
	for beg in mapP[chr]:
	    other=beg+flen
	    s1=0
	    for x in range(other,other+10):
	    	if x in mapM[chr]:
	    	    s1+=mapM[chr][x]
	    if s1>0:
	    	a.append(mapP[chr][beg])
	    	b.append(s1)
	for beg in mapM[chr]:
	    other=beg-flen
	    s1=0
	    for x in range(other,other+10):
	    	if x in mapP[chr]:
	    	    s1+=mapP[chr][x]
	    if s1>0:
	    	b.append(mapM[chr][beg])
	    	a.append(s1)
	corr=np.corrcoef(a,b)[0,1]
	ans+=1.0*cnt[chr]/total*corr
    print '\t'.join([str(flen), str(ans)])


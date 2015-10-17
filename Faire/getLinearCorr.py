import sys,os
import numpy as np

if len(sys.argv)<3:
    print 'Usage: geneExp pathway'
    exit(1)

mapExp={}
mapCnt={}
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[0] in mapExp:
    	mapCnt[temp[0]]+=1
    	mapExp[temp[0]]+=float(temp[1])
    else:
    	mapCnt[temp[0]]=1
    	mapExp[temp[0]]=float(temp[1])
fi.close()

a=[]
b=[]
fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    expLst=[]
    temp=line.split('\t')
    if len(temp)<4 or temp[2] not in mapExp or temp[3] not in mapExp:
    	continue
    a.append(mapExp[temp[2]]/mapCnt[temp[2]])
    b.append(mapExp[temp[3]]/mapCnt[temp[3]])
fi.close()
print len(a), np.corrcoef(a,b)[0][1]

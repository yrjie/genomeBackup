import sys,os
import numpy as np

if len(sys.argv)<3:
    print 'Usage: pos.tab all.tab [mode]'
    exit()

posArr=[]
try:
    mode=int(sys.argv[3])
except:
    mode=0
num=0
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=np.array([float(x) for x in line.split()])
    if num==0:
    	num=len(temp)
    temp/=(sum(temp)+0.1)
    posArr.append(temp)
fi.close()
meanP=np.mean(posArr,axis=0)
if mode>0:
    meanP=np.ones(num)/num
#if mode>0:
#    print meanP
#    exit()
fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=np.array([float(x) for x in line.split()])
    temp/=(sum(temp)+0.1)
    ret=np.sqrt(sum((temp-meanP)**2))
    print ret
fi.close()

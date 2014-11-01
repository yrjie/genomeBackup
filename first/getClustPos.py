import sys,os
import scipy

if len(sys.argv)<3:
    print 'Usage: python getClustPos.py clustFile PosFile'
    exit(1)

class Clust:
    score=0
    allM=[]
    plus=0
    minus=0
    def __init__(self,score):
    	self.score=score

mapS={}
cnt=0
ind=0
allC=[]
fiC=open(sys.argv[1],'r')
fiP=open(sys.argv[2],'r')
for line1 in fiC:
    line2=fiP.readline()
    line1=line1.strip()
    line2=line2.strip()
    if len(line1)<1:
    	continue
    temp1=line1.split('\t')
    temp2=line2.split('\t')
    score=float(temp1[3])
    pos=int(temp2[1])-1000
    if temp2[3]=='-':
    	pos*=-1
    if score not in mapS:
    	mapS[score]=cnt
    	ind=cnt
    	cnt+=1
    	tempC=Clust(score)
    	tempC.allM=[]
    	allC.append(tempC)
    else:
    	ind=mapS[score]
    allC[ind].allM.append(pos)
    if temp2[3]=='+':
    	allC[ind].plus+=1
    else:
    	allC[ind].minus+=1
fiC.close()
fiP.close()
allC.sort(lambda x,y: cmp(x.score,y.score))
for c1 in allC:
    print str(c1.score)+'\t'+str(len(c1.allM))+'\t'+str(scipy.mean(c1.allM))+'\t'+str(scipy.std(c1.allM))+'\t'+str(c1.plus)+'\t'+str(c1.minus)

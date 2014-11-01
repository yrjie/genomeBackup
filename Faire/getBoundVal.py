import sys,os

if len(sys.argv)<3:
    print 'Usage: file threForCTCF'
    exit(1)

class TF:
    name=''
    pct=0
    cov=0
    def __init__(self, name, pct, cov):
    	self.name=name
    	self.pct=pct
    	self.cov=cov

    def __lt__(self, other):
    	return -self.pct<-other.pct

allPct={}
th=float(sys.argv[2])
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    allPct[temp[0]]=[float(temp[i]) for i in range(1,len(temp)) if i%2==0]
fi.close()
cut=0
for x in allPct:
    if 'CTCF' not in x:
    	continue
    now=0
    totalC=sum(allPct[x])
    for i in range(len(allPct[x])):
    	now+=allPct[x][i]
	if now/totalC>th:
	    cut=i
	    break
    break
print cut+1
lstTF=[]
for x in allPct:
    try:
    	pct=sum(allPct[x][0:cut+1])/sum(allPct[x])
    except:
    	pct=0
    lstTF.append(TF(x,pct,sum(allPct[x])))
lstTF.sort()
for x in lstTF:
    print '\t'.join([x.name,str(x.pct),str(x.cov)])

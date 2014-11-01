import sys,os

if len(sys.argv)<3:
    print 'Usage: XYdata binNum'
    exit(1)

class Entry:
    x=0
    y=0
    def __init__(self, x, y):
    	self.x=x
    	self.y=y

binN=int(sys.argv[2])
allE=[]

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    allE.append(Entry(float(temp[0]),float(temp[1])))
fi.close()

minX=min([e.x for e in allE])
maxX=max([e.x for e in allE])

binSize=(maxX-minX)/binN
allBin={}

for e in allE:
    bin=int((e.x-minX)/binSize)*binSize+minX
    if bin not in allBin:
    	allBin[bin]=[]
    else:
    	allBin[bin].append(e.y)

allBinKey=allBin.keys()
allBinKey.sort()
for bin in allBinKey:
    print str(bin)+"\t"+"\t".join([str(y) for y in allBin[bin]])

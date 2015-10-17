import sys,os

if len(sys.argv)<3:
    print 'Usage: posTag negTag'
    exit(1)

maxL=1200

posN=[0]*maxL
negN=[0]*maxL

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    l1=int(temp[-2])
    if l1<0:
    	l1=-l1
    if l1>=maxL:
    	continue
    posN[l1]+=1
fi.close()

fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    l1=int(temp[-2])
    if l1<0:
    	l1=-l1
    if l1>=maxL:
    	continue
    negN[l1]+=1
fi.close()

posR=[1.0*x/sum(posN) for x in posN]
negR=[1.0*x/sum(negN) for x in negN]

for i in range(maxL):
    print '\t'.join([str(x) for x in [i,posR[i]-negR[i]]])

import sys,os
import numpy as np

# the 4th field should be the identifier

if len(sys.argv)<2:
    print 'Usage: overlap'
    exit(1)

fi=open(sys.argv[1])
now=''
short=120
binS=50
offset=20
binN=offset*2+1
num=[1]*binN
peakL=0
plusPos=[]
minusPos=[]
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[3]==now:
    	#leng=int(temp[-2])
    	leng=abs(int(temp[-2]))
    	idx=leng/binS
	if idx<-offset:
	    idx=-offset
	elif idx>offset:
	    idx=offset
	idx+=offset
	num[idx]+=1
	if temp[-1]=='+':
	    plusPos.append(int(temp[12]))
	else:
	    minusPos.append(int(temp[13]))
    else:
	if now!='':
	    #print '\t'.join([str(sum(num))]+[str(1.0*x/sum(num)) for x in num])
	    #print '\t'.join([str(1.0*sum(num)/peakL),str(peakL)]+[str(x) for x in num])
	    if len(plusPos)>0:
	    	plusStd=np.std(plusPos)
	    else:
	    	plusStd=-1
	    if len(minusPos)>0:
	    	minusStd=np.std(minusPos)
	    else:
	    	minusStd=-1
	    print '\t'.join([str(1.0*sum(num)),str(peakL), str(plusStd), str(minusStd)])
	    plusPos=[]
	    minusPos=[]
	    #print '\t'.join([str(x) for x in num+[sum(num)]])
	peakL=int(temp[2])-int(temp[1])
    	now=temp[3]
	num=[1]*binN
#print '\t'.join([str(sum(num))]+[str(1.0*x/sum(num)) for x in num])
#print '\t'.join([str(sum(num)),str(peakL)]+[str(x) for x in num])
if len(plusPos)>0:
    plusStd=np.std(plusPos)
else:
    plusStd=-1
if len(minusPos)>0:
    minusStd=np.std(minusPos)
else:
    minusStd=-1
#print '\t'.join([str(1.0*sum(num)),str(peakL), str(plusStd), str(minusStd)])
print '\t'.join([str(1.0*sum(num)),str(peakL), str(plusStd), str(minusStd)])
#print '\t'.join([str(x) for x in num+[sum(num)]])
fi.close()

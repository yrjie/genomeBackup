import sys,os
import scipy

if len(sys.argv)<2:
    print 'Usage: sigfile'
    exit(1)

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    num=[float(x) for x in temp]
    ind1=num.index(max(num))
    a=max(0,ind1-20)
    b=ind1+20
    ind2=num.index(min(num[a:b]))
    #print str(num[ind1]-num[ind2])
    #print str(ind2-ind1)+"\t"+str(num[ind1]-num[ind2])
    #print scipy.std(num)/(scipy.median(num)+0.1)
    print scipy.std(num)/(scipy.mean(num)+0.1)
    #print scipy.mean(num)+0.1
    #print scipy.log(sum(num)+0.1)
fi.close()

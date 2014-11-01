import sys,os
import scipy.stats as ss

if len(sys.argv)<2:
    print 'Usage: dat1 [dat2]'
    exit(1)

a=[]
b=[]

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    a.append(float(temp[0]))
    if len(temp)>1:
    	b.append(float(temp[1]))
fi.close()

if len(sys.argv)>2:
    fi=open(sys.argv[2])
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split()
	b.append(float(temp[0]))
    fi.close()

print ss.spearmanr(a,b)

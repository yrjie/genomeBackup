import sys,os

if len(sys.argv)<2:
    print 'Usage: python checkPattern.py patternCount.file [threshold]'
    exit(1)

th=1
if len(sys.argv)>=3:
    th=float(sys.argv[2])

fi=open(sys.argv[1])
cnt=[0 for i in range(4)]
all=0
numP=0
ok=0
fi.readline()
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    for j in range(4):
    	cnt[j]+=int(temp[j])
    	all+=int(temp[j])
    if int(temp[1])>th and int(temp[2])>th:
    	ok+=1
    numP+=1
fi.close()
for i in range(4):
    print str(cnt[i])+'\t'+str(1.0*cnt[i]/all)
print 1.0*all/4/numP
print 1.0*ok/numP

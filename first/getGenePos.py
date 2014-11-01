import sys,os

if len(sys.argv)<3:
    print 'Usage: gene.lst gene.bed'
    exit(1)

mapT=set()
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    mapT.add(line)
fi.close()
fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[4] in mapT:
    	print line
fi.close()

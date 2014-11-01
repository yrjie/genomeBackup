import sys,os

if len(sys.argv)<2:
    print 'Usage: python intersectGene.py infile'
    exit(1)

fi=open(sys.argv[1],'r')
a1=set()
a2=set()
b=a1
for line in fi:
    line=line.strip()
    if len(line)<1:
    	b=a2
    	continue
    b.add(line)
fi.close()
b=a1&a2
print str(len(a1))+"\t"+str(len(a2))
print str(len(b))+"\t"+str(b)

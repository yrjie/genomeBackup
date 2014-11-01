import sys,os

if len(sys.argv)<3:
    print 'Usage: smpID mapping'
    exit(1)

mapG={}
fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    mapG[temp[1]]=temp[0]
fi.close()

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    print mapG[line]
fi.close()

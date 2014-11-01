import sys,os

if len(sys.argv)<4:
    print 'Usage: chrLen chrCnt gapfile'
    exit(1)

chrLen={}
chrCnt={}

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    chrLen[temp[0]]=temp[1]
fi.close()

fi=open(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    chrCnt[temp[1]]=temp[0]
fi.close()

fi=open(sys.argv[3])
nowChr=''
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    if temp[0] != nowChr:
    	nowChr=temp[0]
    	if temp[0] in chrLen and temp[0] in chrCnt:
    	    print '\t'.join([nowChr, chrLen[nowChr], chrCnt[nowChr]])
    print line
fi.close()


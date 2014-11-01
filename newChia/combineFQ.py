import string
import sys,os

if len(sys.argv)<3:
    print 'Usage: seqFile lenOfSubStr'
    exit(1)

def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

def comple(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    cseq = dna.translate(complements)
    return cseq

lenSub=int(sys.argv[2])
allRead={}
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split()
    if len(temp)>2:
    	if temp[0]!=temp[2]:
    	    raise Exception('Non-equal ID: '+temp[0])
    	id=temp[0]
    	continue
    a=temp[0]
    b=rc(temp[1])
    ainb=a[-lenSub:] in b
    bina=b[-lenSub:] in a #or b[0:lenSub] in a
    #print a,b
    if ainb:
    	print '>'+id
    	print a+b[b.index(a[-lenSub:])+lenSub:]
    elif bina:
    	print '>'+id
    	ret=b+a[a.index(b[-lenSub:])+lenSub:]
	if a[0:10] in ret:
	    head=ret.index(a[0:10])
	    tail=ret.index(b[-lenSub:])+lenSub
	    print ret[head:tail]
	else:
	    print ret
fi.close()


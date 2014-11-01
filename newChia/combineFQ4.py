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
while 1:
    l1=fi.readline().strip()
    if len(l1)<1:
    	break
    l2=fi.readline().strip()
    l3=fi.readline().strip()
    l4=fi.readline().strip()
    t1=l1.split()
    t2=l2.split()
    t3=l3.split()
    t4=l4.split()
    if t1[0]!=t1[2]:
	raise Exception('Non-equal ID: '+t1[0])
    id=t1[0]
    a=t2[0]
    aq=t4[0]
    b=rc(t2[1])
    bq=t4[1][::-1]
    ainb=a[-lenSub:] in b
    bina=b[-lenSub:] in a #or b[0:lenSub] in a
    #print a,b
    if ainb:
    	print id
    	print a+b[b.index(a[-lenSub:])+lenSub:]
    	print '+'
	print aq+bq[b.index(a[-lenSub:])+lenSub:]
    elif bina:
    	print id
    	ret=b+a[a.index(b[-lenSub:])+lenSub:]
	retq=bq+aq[a.index(b[-lenSub:])+lenSub:]
	if a[0:10] in ret:
	    head=ret.index(a[0:10])
	    tail=ret.index(b[-lenSub:])+lenSub
	    ret=ret[head:tail]
	    retq=retq[head:tail]
	print ret
	print '+'
	print retq
fi.close()


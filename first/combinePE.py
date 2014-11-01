import string
import sys,os

if len(sys.argv)<3:
    print 'Usage: samfile lenOfSubStr'
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
    temp=line.split('\t')
    flag=int(temp[1])
    strand=flag&32
    num=flag&64
    if len(temp[9])<100:
    	continue
    if strand==0:
    	seq=rc(temp[9])
    else:
    	seq=temp[9]
    now=temp[0]+"_"+str(num)
    if num==0:
    	seq=rc(seq)
    if now not in allRead:
    	allRead[now]=seq
fi.close()

ret={}
for x in allRead:
    temp=x.split('_')
    name=temp[0]
    num=int(temp[1])
    if name not in ret:
#    	ret[name]='\t'.join([name,allRead[x]])
	ret[name]=allRead[x]
    else:
    	ret[name]+='\t'+allRead[x]
    	a=ret[name].split('\t')[0]
    	b=ret[name].split('\t')[1]
	ainb=a[-lenSub:] in b #or a[0:lenSub] in b
	bina=b[-lenSub:] in a #or b[0:lenSub] in a
	#print a,b
	if ainb:
	    print a+b[b.index(a[-lenSub:])+lenSub:]
	elif bina:
	    print b+a[a.index(b[-lenSub:])+lenSub:]
    	#print '\t'.join([name, a,b,str(ainb),str(bina)])

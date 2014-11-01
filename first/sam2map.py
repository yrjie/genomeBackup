import sys,os

if len(sys.argv)<2:
    print 'Usage: samfile \n The samfile should be sorted'
    exit(1)

class Tag:
    name=''
    chr=''
    pos=0
    strand=''
    seq=''
    slen=0
    isMul=False
    def __init__(self, name='', chr='', pos=0, strand='', seq='', isMul=False):
    	self.slen=len(seq)
    	self.name=name
    	self.chr=chr
    	self.pos=pos
    	self.strand=strand
    	self.seq=seq
    	self.isMul=isMul
#	if self.strand=='-':
#	    self.pos+=(self.slen-1)

    def __eq__(self, other):
    	return self.chr==other.chr and self.pos==other.pos and self.slen==other.slen

    def __ne__(self, other):
    	return not (self==other)

def printTag(tag, type):
    pos=tag.pos
    if tag.strand=='-':
    	pos+=tag.slen-1
    print '\t'.join([type, tag.chr, tag.strand, str(pos), '0', str(tag.slen)])
    if tag.isMul:
    	print '\t'.join([type, tag.chr, tag.strand, str(pos), '0', str(tag.slen)])

allTag=[]
tagMap={}
fi=open(sys.argv[1])
cnt=0
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    name=temp[0]
    chr=temp[2]
    pos=int(temp[3])
    seq=temp[9]
    isMul=int(temp[4])==0
    if int(temp[1])&16>0:
    	strand='-'
    else:
    	strand='+'
    allTag.append(Tag(name, chr,pos,strand,seq, isMul))
    tagMap[temp[0]]=cnt
    cnt+=1
fi.close()

pre1=Tag()
pre2=Tag()
cnt=0
for idx, tag1 in enumerate(allTag):
    if tag1.name[-1]=='2':
    	tmpname=list(tag1.name)
    	tmpname[-1]='1'
    	name1=''.join(tmpname)
	if name1 not in tagMap:
	    print '\t'.join(['>'+tag1.name[:-2]+'_COUNT:1',tag1.seq, tag1.seq])
	    printTag(tag1, 'T')
	continue

    tmpname=list(tag1.name)
    tmpname[-1]='2'
    name2=''.join(tmpname)
    if name2 not in tagMap:
	print '\t'.join(['>'+tag1.name[:-2]+'_COUNT:1',tag1.seq, tag1.seq])
	printTag(tag1, 'H')
    	continue
    tag2=allTag[tagMap[name2]]
    if pre1!=tag1 or pre2!=tag2:
#	if cnt>1:
#	    print 'cnt',pre1.name
	if pre1.name!='':
	    print '\t'.join(['>'+pre1.name[:-2]+'_COUNT:'+str(cnt),pre1.seq, pre2.seq])
	    printTag(pre1, 'H')
	    printTag(pre2, 'T')
#	    pos=pre1.pos
#	    if pre1.strand=='-':
#		pos+=(pre1.slen-1)
#	    print '\t'.join(['H', pre1.chr, pre1.strand, str(pos), '0', str(pre1.slen)])
#	    
#	    pos=pre2.pos
#	    if pre2.strand=='-':
#		pos+=(pre2.slen-1)
#	    print '\t'.join(['T', pre2.chr, pre2.strand, str(pos), '0', str(pre2.slen)])
    	pre1=tag1
    	pre2=tag2
    	cnt=1
    else:
    	cnt+=1

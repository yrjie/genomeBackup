import sys,os
import resource

if len(sys.argv)<2:
    print 'Usage: samfile\nThe file should be sorted by ID and all ID are assumed to have ID_1 and ID_2'
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
    
def temp2tag(temp):
    name=temp[0]
    chr=temp[2]
    pos=int(temp[3])
    seq=temp[9]
    isMul=int(temp[4])==0
    if int(temp[1])&16>0:
        strand='-'
    else:
        strand='+'
    return Tag(name, chr, pos, strand, seq, isMul)

def printTag(tag, type):
    if tag.chr=='*':
        return
    pos=tag.pos
    if tag.strand=='-':
    	pos+=tag.slen-1
    print '\t'.join([type, tag.chr, tag.strand, str(pos), '0', str(tag.slen)])
    if tag.isMul:
    	print '\t'.join([type, tag.chr, tag.strand, str(pos), '0', str(tag.slen)])


mapTag={}
fi=open(sys.argv[1])
cnt=0
while True:
    line1=fi.readline()
    line1=line1.strip()
    if len(line1)<1:
    	break
    temp1=line1.split('\t')
    id1=temp1[0][:-2]
    line2=fi.readline()
    temp2=line2.split('\t')
    id2=temp2[0][:-2]
    if id1!=id2:
        raise Exception(id1)
    chr1=temp1[2]
    pos1=temp1[3]
    slen1=str(len(temp1[9]))
    chr2=temp2[2]
    pos2=temp2[3]
    slen2=str(len(temp2[9]))
    key='_'.join([chr1,pos1,slen1,chr2,pos2,slen2])
    if key not in mapTag:
    	mapTag[key]=[0, False]
    mapTag[key][0]+=1
    cnt+=1
    if cnt%10000000==0:
        print >>sys.stderr, '\t'.join([str(cnt),str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)])
fi.close()

fi=open(sys.argv[1])
while True:
    line1=fi.readline()
    line1=line1.strip()
    if len(line1)<1:
        break
    temp1=line1.split('\t')
    id1=temp1[0][:-2]
    line2=fi.readline()
    temp2=line2.split('\t')
    id2=temp2[0][:-2]
    assert id1==id2
    chr1=temp1[2]
    pos1=temp1[3]
    slen1=str(len(temp1[9]))
    chr2=temp2[2]
    pos2=temp2[3]
    slen2=str(len(temp2[9]))
    key='_'.join([chr1,pos1,slen1,chr2,pos2,slen2])
    if not mapTag[key][1]:
        mapTag[key][1]=True
        tag1=temp2tag(temp1)
        tag2=temp2tag(temp2)
        print '\t'.join(['>'+tag1.name[:-2]+'_COUNT:'+str(mapTag[key][0]),tag1.seq, tag2.seq])
        printTag(tag1, 'H')
        printTag(tag2, 'T')
fi.close()


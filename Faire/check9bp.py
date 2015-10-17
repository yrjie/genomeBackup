import sys,os

if len(sys.argv)<2:
    print 'Usage: sortedTag [printTag]'
    exit(1)

class Tag:
    chr=''
    start=0
    end=0
    fragl=0
    strand='+'
    def __init__(self, chr, start, end, fragl, strand):
    	self.chr=chr
    	self.start=start
    	self.end=end
	if fragl<0:
	    fragl=-fragl
    	self.fragl=fragl
    	self.strand=strand

def overlap(tag1, tag2):
    if tag1.chr!=tag2.chr or tag1.end<=tag2.start:
    	return 0
    return tag1.end-tag2.start

allTag=[]
numTag=0

ptag=0
if len(sys.argv)>2:
    ptag=int(sys.argv[2])

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    allTag.append(Tag(temp[0],int(temp[1]),int(temp[2]),int(temp[4]),temp[5]))
    numTag+=1
fi.close()
#print numTag

i=0
while i<numTag:
    if allTag[i].strand=='+':
    	i+=1
    	continue
    j=i+1
    while 1:
	if j>=numTag:
	    break
    	obp=overlap(allTag[i], allTag[j])
	if obp==0:
	    break
	if allTag[j].strand=='+':
	    if not ptag:
	    	print obp
	    elif obp==9:
	    	#print '\t'.join([allTag[i].chr,str(allTag[j].start),str(allTag[i].end)])
		if allTag[i].fragl<600 and allTag[j].fragl<600:
		    print '\t'.join([str(allTag[i].fragl),str(allTag[j].fragl)])
    	j+=1
    i=j


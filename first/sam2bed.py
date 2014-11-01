import sys,os

if len(sys.argv)<2:
    print 'Usage: python sam2bed.py samfile [outfile]\nIf outfile is not provided, \'out.bed\' will be used'
    exit(1)

outfile='out.bed'
if len(sys.argv)>=3:
    outfile=sys.argv[2]

def getStrand(samFlag):
    strand='+'
    if samFlag&(0x10):
    	strand='-'
    return strand

fi=open(sys.argv[1],'r')
fo=open(outfile,'w')
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    strand=getStrand(int(temp[1]))
    end=int(temp[3])+len(temp[9])+1
    fo.write(temp[2]+'\t'+temp[3]+'\t'+str(end)+'\t'+strand+'\n')
fi.close()
fo.close()

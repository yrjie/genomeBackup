import sys,os

if len(sys.argv)<2:
    print 'Usage: python getTss.py geneFile'
    exit(1)

fi=open(sys.argv[1],'r')
out=range(0,4)
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[5]=='+':
    	tss=int(temp[1])
    else:
    	tss=int(temp[2])
    out[0]=temp[0]
    out[1]=str(tss)
    out[2]=str(tss)
    out[3]=temp[3]+';'+temp[4]
    print '\t'.join(out)
fi.close()

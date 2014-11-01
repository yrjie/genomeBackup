import sys,os

if len(sys.argv)<3:
    print 'Usage: python binChr.py genomefile binsize'
    exit(1)

fi=open(sys.argv[1])
bin=int(sys.argv[2])
for line in fi:
    line=line.strip()
    if len(line)<1 or '_' in line:
    	continue
    temp=line.split('\t')
    start=1
    chrlen=int(temp[1])
    while start<chrlen:
    	print temp[0]+'\t'+str(start)+'\t'+str(min(start+bin,chrlen))
    	start+=bin
fi.close()

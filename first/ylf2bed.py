import sys,os

if len(sys.argv)<3:
    print 'Usage: python ylf2bed.py ylfFile chromInfo'
    exit(1)

chrlen={}
readlen=36
fi=open(sys.argv[2],'r')
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    chrlen[temp[0][3:]]=int(temp[1])
fi.close()

fi=open(sys.argv[1],'r')
now=''
cnt=0
strand=''
for line in fi:
    line=line.strip()
    if line[0]=='>':
	now=line[1:]
	continue
    try:
    	start=int(line)
    except:
    	continue
    strand='+'
    if start<0:
    	start=-start#+chrlen[now]
    	strand='-'
    end=start+36
    cnt+=1
    print 'chr'+now+'\t'+str(start)+'\t'+str(end)+'\tbowtie'+str(cnt)+'\t1\t'+strand
fi.close()

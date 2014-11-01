import sys,os

class Peak:
    chr='chr1'
    start='0'
    end='0'
    score=0
    def __init__(self,ch,a,b,c):
    	self.chr=ch
    	self.start=a
    	self.end=b
    	self.score=c

def filter(infile, top):
    fi=open(infile,'r')
    fo=open('filtered.narrowPeak','w')
    all=[]
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	p1=int(temp[9])
	if p1==-1:
	    p1=(int(temp[1])+int(temp[2]))/2
	else:
	    p1=int(temp[1])+p1
	a=p1-200
	b=p1+200
	t1=Peak(temp[0],str(a),str(b),int(temp[4]))
	all.append(t1)
    all.sort(lambda x,y:-cmp(x.score,y.score))
    #print all[0].score
    cnt=0
    for pk in all:
    	fo.write(pk.chr+'\t'+pk.start+'\t'+pk.end+'\t.\t'+str(pk.score)+'\t.\n')
    	cnt+=1
	if cnt>=top:
	    break
    fi.close()
    fo.close()

def getCov(inBW, inCH, th):
    cov=0
    cmd='/home/ruijie/bin/bigWigSummaryBatch '+inBW+' '+inCH+' 1 >temp'
    os.system(cmd)
    ftemp=open('temp','r')
    for line in ftemp:
    	line=line.strip()
	if len(line)<1:
	    continue
	ans=float(line)
	if ans>th:
	    cov+=1
    ftemp.close()
    return cov


if len(sys.argv)<4:
    print 'Usage: python dnaseCovBW.py DNase.bigWig chipFileList threshold [outfile]\nIf outfile is not provided, <prefix>.cov will be used'
    exit(1)

chipPath='/home/browser/BASIC/import-data/narrowPeak/'
temp=sys.argv[1].split('/')
outfile='coverage/'+temp[len(temp)-1].split('.')[0]+'.cov'
th=float(sys.argv[3])
if len(sys.argv)>4:
    outfile=sys.argv[4]
print outfile

fi=open(sys.argv[2],'r')
fo=open(outfile,'w')

for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    infile=chipPath+line
#    pNum=len(open(infile,'r').readlines())
#    cov=getCov(sys.argv[1],infile,th)
    filter(infile,10000)
    pNum=len(open('filtered.narrowPeak','r').readlines())
    cov=getCov(sys.argv[1],'filtered.narrowPeak',th)
    print cov
#    cmd='bedtools intersect -a '+sys.argv[1]+' -b '+line+' >temp'
#    print cmd
#    os.system(cmd)
#    cov=len(open('temp','r').readlines())
#    cmd='bedtools shuffle -i '+infile+' -g ~/genome/human.hg19.genome >shuffled'
    cmd='bedtools shuffle -i filtered.narrowPeak -g ~/genome/human.hg19.genome >shuffled'
#    print cmd
    os.system(cmd)
    FP=getCov(sys.argv[1],'shuffled',th)
    print FP
#    cmd='bedtools intersect -a '+sys.argv[1]+' -b shuffled >tempS'
#    print cmd
#    os.system(cmd)
#    FP=len(open('tempS','r').readlines())
    line=line.split('.')[0]
#    print 1.0*FP/pNum
    fo.write(str(pNum)+'\t'+str(1.0*cov/pNum)+'\t'+str(1.0*FP/pNum)+'\t'+line+'\n')

fi.close()
fo.close()

import sys,os
import scipy

if len(sys.argv)<3:
    print 'Usage: python getClustBed.py clustFile prefix'
    exit(1)

class Region:
    chr=''
    start=''
    end=''
    score=0
    def __init__(self,chr,start,end,score):
    	self.chr=chr
    	self.start=start
    	self.end=end
    	self.score=score

prefix=sys.argv[2]+'_c'
allC=[]
fiC=open(sys.argv[1],'r')
for line in fiC:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    allC.append(Region(temp[0],temp[1],temp[2],float(temp[3])))
fiC.close()
allC.sort(lambda x,y: cmp(x.score,y.score))
now=0
cnt=0
for c1 in allC:
    if now!=c1.score:
    	cnt+=1
    	now=c1.score
	try:
	    fo.close()
	except:
	    1
    	fo=open(prefix+str(cnt)+'.bed','w')
    fo.write(c1.chr+'\t'+c1.start+'\t'+c1.end+'\t'+str(c1.score)+'\n')
fo.close()

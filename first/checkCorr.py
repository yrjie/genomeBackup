import sys,os
import numpy as np

if len(sys.argv)<4:
    print 'Usage: python checkCorr.py DNaseFile peakList threshold [window] [top]\nIf window is not provided, 200 will be used\nIf top is not provided, 10000 will be used'
    exit(1)

win=200
if len(sys.argv)>=5:
    win=int(sys.argv[4])

top=10000
if len(sys.argv)>=6:
    top=int(sys.argv[5])

class Peak:
    chr='chr1'
    start='0'
    end='0'
    val=0
    def __init__(self,ch,a,b,c):
    	self.chr=ch
    	self.start=a
    	self.end=b
    	self.val=c

def correlation(t1,t2,th):
    data1=np.array(map(float,open(t1,'r').readlines()))
    data2=np.array(map(float,open(t2,'r').readlines()))
    idx=np.nonzero(data1>th)
    cover=1.0*len(idx[0])/len(data1)
    data1=data1[idx]
    data2=data2[idx]
    m1=data1.mean()
    m2=data2.mean()
    std1=data1.std()
    std2=data2.std()
    cov=((data1-m1)*(data2-m2)).mean()
#    print m1
#    print m2
#    print std1
#    print std2
#    print cov
    corr=cov/(std1*std2)
    print str(corr)+'\t'+str(cover)
    return [corr,cover]

def filter(infile, win, top):
    fi=open(infile,'r')
    fo=open('filtered.peak','w')
    all=[]
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	p1=-1
	if len(temp)>=10:
	    p1=int(temp[9])
	if p1==-1:
	    p1=(int(temp[1])+int(temp[2]))/2
	else:
	    p1=int(temp[1])+p1
	a=p1-win
	b=p1+win
#	a=int(temp[1])
#	b=int(temp[2])
	t1=Peak(temp[0],str(a),str(b),float(temp[6]))
	all.append(t1)
    all.sort(lambda x,y:-cmp(x.val,y.val))
    #print all[0].score
    cnt=0
    for pk in all:
    	fo.write(pk.chr+'\t'+pk.start+'\t'+pk.end+'\t.\t.\t.\t'+str(pk.val)+'\n')
    	cnt+=1
	if cnt>=top:
	    break
    fi.close()
    fo.close()

chipPath='/home/browser/BASIC/import-data/'
fiL=open(sys.argv[2],'r')
th=float(sys.argv[3])
outfile='correlation/'+os.path.basename(sys.argv[1]).split('.')[0]+'.corr'
print outfile
foCorr=open(outfile,'w')
for line in fiL:
    line=line.strip()
    if len(line)<1 or line[0]=='#':
    	continue
    if '\t' in line:
    	temp=line.split('\t')
    else:
    	temp=line.split('.')
    inCH=chipPath+temp[1]+'/'+temp[0]+'.'+temp[1]
    filter(inCH,win,top)
    inCH='filtered.peak'
    cmd='/home/ruijie/bin/bigWigSummaryBatch '+sys.argv[1]+' '+inCH+' 1 >val_DNase'
    os.system(cmd)
    fiTF=open(inCH,'r')
    foTF=open('val_TF','w')
    for item in fiTF:
    	temp=item.strip().split('\t')
	if len(temp)<7:
	    continue
	foTF.write(temp[6]+'\n')
    foTF.close()
    ans=correlation('val_DNase','val_TF',th)
    foCorr.write(str(ans[0])+'\t'+str(ans[1])+'\t'+line.split('.')[0]+'\n')
fiL.close()
foCorr.close()


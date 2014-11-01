import os,sys

if len(sys.argv)<4:
    print 'Usage: python train.py TF_name tfbs.lst dnase.lst [threshold_DNase]'
    exit(1)

thD=1.7
if len(sys.argv)>=5:
    thD=float(sys.argv[4])

class TFBS:
    file=''
    cell=''
    def __init__(self,f1, c1):
    	self.file=f1
    	self.cell=c1

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

def read_tfbs(tf_name,tfbs_file,cellS_tf):
    dir0='./data/'
    ret=[]
    fi=open(tfbs_file,'r')
    for line in fi:
    	temp=line.split('\t')
    	temp[0]=dir0+temp[0]+'.narrowPeak'
	if tf_name+' ' in temp[1] or tf_name+'\n' in temp[1]:
	    cell=temp[1].split(' ')[0]
	    if '+' in cell:
	    	cell=cell.split('+')[0]
	    cellS_tf.add(cell)
	    obj=TFBS(temp[0],cell)
	    ret.append(obj)
    fi.close()
#    print len(ret)
    return ret

def read_dnase(infile):
    fi=open(infile,'r')
    ret={}
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	cell=temp[1][0:temp[1].find(' OS')]
	if cell in ret:
	    sys.stderr.write(cell+'\n')
	ret[cell]=temp[0]
    fi.close()
    return ret

def union(cmdU):
    cmdU+=' >temp'
    os.system(cmdU)
    os.system('bedtools sort -i temp >tempS')
    os.system('bedtools merge -i tempS -nms >tempM')
    return 'tempM'

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
	if win>0:
	    a=p1-win
	    b=p1+win
	else:
	    a=int(temp[1])
	    b=int(temp[2])
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

def coverage(afile, bfile):
    covFile='tf.cov'
    cmd='bedtools coverage -a '+afile+' -b '+bfile+' >'+covFile
    th=0.5
    os.system(cmd)
    aNum=len(open(afile,'r').readlines())
    bNum=len(open(bfile,'r').readlines())
    yes=0
    fi=open(covFile,'r')
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	cov1=float(temp[len(temp)-1])
	if cov1>th:
	    yes+=1
    fi.close()
    return [1.0*yes/aNum,1-1.0*yes/bNum]

def predict(uniFile, dnaseFile, th):
    dnaseFile='/home/ruijie/gbdb/hg19/bbi/'+dnaseFile+'.bigWig'
    sigFile='predicted.sig'
    peakFile='predicted.peak'
    cmd='bigWigSummaryBatch '+dnaseFile+' '+uniFile+' 1 >'+sigFile
    os.system(cmd)
    fi=open(sigFile,'r')
    fiU=open(uniFile,'r')
    fo=open(peakFile,'w')
    for line in fi:
    	pk=fiU.readline()
    	line=line.strip()
	if len(line)<1:
	    continue
	if float(line)>th:
	    fo.write(pk+'\n')
    fi.close()
    fiU.close()
    fo.close()
    return peakFile

def work(tfbs,cellS_tf,dnaseDict):
    name=[]
    ans=[]
    fpr=[]
    bNum=[]
    uniFile=''
    if len(cellS_tf)==1:
    	ans.append(0)
    	ld=len(ans)
    	print sys.argv[1]+'\t'+str(ans[0])+'\t'+str(ans[ld/4])+'\t'+str(ans[ld/2])+'\t'+str(ans[3*ld/4])+'\t'+str(ans[ld-1])
    	return
    for cell in cellS_tf:
	if cell not in dnaseDict:
	    sys.stderr.write('not have dnase: '+cell+'\n')
	    continue
    	idx=[]
    	cmdU='cat'
	for i in range(0,len(tfbs)):
	    if tfbs[i].cell!=cell:
	    	cmdU+=' '+tfbs[i].file
	    else:
	    	idx.append(i)
	uniFile=union(cmdU)
	ma=0
	fp1=0
	na1=''
	preFile=predict(uniFile, dnaseDict[cell], thD)
	bNum.append(len(open(preFile,'r').readlines()))
	for i in idx:
	    filter(tfbs[i].file,-1,10000)
	    cov=coverage('filtered.peak', preFile)
	    #cov=coverage(tfbs[i].file, preFile)
	    if cov[0]>ma:
	    	na1=tfbs[i].file
	    	ma=cov[0]
	    	fp1=cov[1]
	ans.append(ma)
	fpr.append(fp1)
	name.append(na1)
    for i in range(len(ans)):
    	print name[i]+'\t'+str(ans[i])+'\t'+str(fpr[i])+'\t'+str(bNum[i])
    ld=len(ans)
    ans.sort()
#    print sys.argv[1]+'.'+str(ld)+'\t'+str(ans[0])+'\t'+str(ans[ld/4])+'\t'+str(ans[ld/2])+'\t'+str(ans[3*ld/4])+'\t'+str(ans[ld-1])

cellS_tf=set()
tfbs=read_tfbs(sys.argv[1],sys.argv[2],cellS_tf)
dnaseDict=read_dnase(sys.argv[3])
work(tfbs,cellS_tf,dnaseDict)

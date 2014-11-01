import os,sys

if len(sys.argv)<3:
    print 'Usage: python checkTFCov.py TF_name tfbs.lst'
    exit(1)

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

def read_tfbs(tf_name,tfbs_file,cellS):
    dir0='/home/browser/BASIC/import-data/narrowPeak/'
    ret=[]
    fi=open(tfbs_file,'r')
    for line in fi:
    	temp=line.split('\t')
    	temp[0]=dir0+temp[0]+'.narrowPeak'
	if tf_name+' ' in temp[1] or tf_name+'\n' in temp[1]:
	    cell=temp[1].split(' ')[0]
	    if '+' in cell:
	    	cell=cell.split('+')[0]
	    cellS.add(cell)
	    obj=TFBS(temp[0],cell)
	    ret.append(obj)
    fi.close()
#    print len(ret)
    return ret

def union(cmdU):
    cmdU+=' >temp'
    os.system(cmdU)
    os.system('bedtools sort -i temp >tempS')
    os.system('bedtools merge -i tempS >tempM')
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
    all=0
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
	all+=1
    fi.close()
    return 1.0*yes/all

def work(tfbs,cellS):
    ans=[]
    uniFile=''
    if len(cellS)==1:
    	ans.append(0)
    	ld=len(ans)
    	print sys.argv[1]+'\t'+str(ans[0])+'\t'+str(ans[ld/4])+'\t'+str(ans[ld/2])+'\t'+str(ans[3*ld/4])+'\t'+str(ans[ld-1])
    	return
    for cell in cellS:
    	idx=[]
    	cmdU='cat'
	for i in range(0,len(tfbs)):
	    if tfbs[i].cell!=cell:
	    	cmdU+=' '+tfbs[i].file
	    else:
	    	idx.append(i)
	uniFile=union(cmdU)
	ma=0
	for i in idx:
	    filter(tfbs[i].file,-1,10000)
	    ma=max(ma,coverage(uniFile,'filtered.peak'))
#	    ma=max(ma,coverage(uniFile,tfbs[i].file))
	ans.append(ma)
    #print ans
    ld=len(ans)
    ans.sort()
    print sys.argv[1]+'.'+str(ld)+'\t'+str(ans[0])+'\t'+str(ans[ld/4])+'\t'+str(ans[ld/2])+'\t'+str(ans[3*ld/4])+'\t'+str(ans[ld-1])

cellS=set()
tfbs=read_tfbs(sys.argv[1],sys.argv[2],cellS)
work(tfbs,cellS)

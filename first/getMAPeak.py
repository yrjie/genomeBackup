import sys,os

if len(sys.argv)<4:
    print 'Usage: python getMAPeak.py MA_result probe_desc microarray_data'
    exit(1)

class Probe:
    name=''
    chr=''
    beg=0
    end=0
    A_val=0
    B_val=0
    def __init__(self, name='', chr='', beg=0, end=0):
    	self.name=name
    	self.chr=chr
    	self.beg=beg
    	self.end=end

mapP={}
allP=[]

def read_probe(infile):
    fi=open(infile,'r')
    fi.readline()
    cnt=0
    ext=2500
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	if '_' in temp[1] or temp[0] in mapP:
	    continue
	if temp[4]=='+':
	    tss=int(temp[2])
	elif temp[4]=='-':
	    tss=int(temp[3])
	beg=tss-ext
	end=tss+ext
	pro1=Probe(temp[0],temp[1],beg,end)
	allP.append(pro1)
	mapP[temp[0]]=cnt
	cnt+=1
    fi.close()

def read_array(infile):
    fi=open(infile,'r')
    fi.readline()
    for line in fi:
	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	if temp[0] not in mapP:
	    continue
	allP[mapP[temp[0]]].A_val=float(temp[1])
	allP[mapP[temp[0]]].B_val=float(temp[2])
    fi.close()

def work(infile):
    fi=open(infile,'r')
    fo=open('out.xls','w')
    fi.readline()
    th_M=0.8
    th_fold=1.5
    fo.write('chr\tstart\tend\tprobe_id\tintensity_A\tintensity_B\n')
    prom=0
    num=0
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	if abs(float(temp[6]))<th_M:
	    continue
	num+=1
	ok=0
	isP=0
	pro1=Probe()
	beg=int(temp[1])
	end=int(temp[2])
	for probe in allP:
	    if probe.chr!=temp[0] or probe.end<beg or probe.beg>end:
	    	continue
	    isP=1
	    if probe.B_val>probe.A_val*th_fold:
	    	ok=1
	    	pro1=probe
	    	break
	if ok:
	    fo.write(temp[0]+'\t'+temp[1]+'\t'+temp[2]+'\t'+pro1.name+'\t'+str(pro1.A_val)+'\t'+str(pro1.B_val)+'\n')
	prom+=isP
    print num
    print prom
    fi.close()
    fo.close()

read_probe(sys.argv[2])
read_array(sys.argv[3])
work(sys.argv[1])

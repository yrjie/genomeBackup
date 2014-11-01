import sys,os
import math

if len(sys.argv)<4:
    print 'Usage: python getMAPlot.py MA_result probe_desc microarray_data'
    exit(1)

class Probe:
    name=''
    chr=''
    beg=0
    end=0
    A_val=1
    B_val=1
    gene=''
    def __init__(self, name='', chr='', beg=0, end=0, gene=''):
    	self.name=name
    	self.chr=chr
    	self.beg=beg
    	self.end=end
    	self.gene=gene

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
	    beg=tss-ext
	    end=tss
	elif temp[4]=='-':
	    tss=int(temp[3])
	    beg=tss
	    end=tss+ext
#	beg=tss-ext
#	end=tss+ext
	pro1=Probe(temp[0],temp[1],beg,end,temp[5])
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
    fo=open('MAPlot.dat','w')
    foR=open('JARID_region.xls','w')
#    fi.readline()
    th_M=0.8
    th_fold=1.5
#    fo.write('chr\tstart\tend\tprobe_id\tintensity_A\tintensity_B\n')
    prom=0
    num=0
    TP=0
    FP=0
    foR.write('chr\tstart\tend\tconfidence value\trelated gene\tlog(intensity_A)\tlog(intensity_B)\tlog(fold change)\n')
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	temp=line.split('\t')
	try:
	    if len(temp)>6 and abs(float(temp[6]))<th_M:
	    	continue
	except:
	    continue
	num+=1
	ok=0
	isP=0
	pro1=Probe()
	try:
	    beg=int(temp[1])
	    end=int(temp[2])
	except ValueError:
	    continue
	for probe in allP:
	    if probe.chr!=temp[0] or probe.end<beg or probe.beg>end:
	    	continue
	    isP=1
	    pro1=probe
	    break
	    #if probe.B_val>probe.A_val*th_fold:
	    #	ok=1
	    #	pro1=probe
	    #	break
	if isP:
	    if 'rma' in sys.argv[3]:
	    	res=pro1.B_val-pro1.A_val
	    else:
	    	res=math.log(pro1.B_val/pro1.A_val)
	    if len(temp)>6:
	    	x=temp[6]
	    else:
	    	x=temp[4]
	    fo.write(x+'\t'+str(res)+'\n')
	    if float(x)>0:
	    	if res>0:
	    	    TP+=1
	    	    foR.write(pro1.chr+'\t'+str(beg)+'\t'+str(end)+'\t'+x+'\t'+pro1.gene+'\t'+str(pro1.A_val)+'\t'+str(pro1.B_val)+'\t'+str(res)+'\n')
		elif res<0:
		    FP+=1
	prom+=isP
    print num
    print 'at promoter: '+str(prom)
    print 'TP: '+str(TP)
    print 'FP: '+str(FP)
    fi.close()
    fo.close()
    foR.close()

read_probe(sys.argv[2])
read_array(sys.argv[3])
work(sys.argv[1])

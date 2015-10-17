import sys,os
import numpy as np
import statsmodels.api as sm
import operator
import copy
import scipy.stats as ss

if len(sys.argv)<3:
    print 'Usage: allPeakTag outfile'
    print 'The 2th ending column should be the fragment length'
    exit(1)

infile=sys.argv[1]
outfile=sys.argv[2]

binS=50
offset=500
binN=2*offset/binS

debug=0
if len(sys.argv)>4:
    debug=int(sys.argv[4])

expC=3
atacC=4
eps=1e-16
allPeak=[]
minB=85
maxB=164

beta=1 # pseudo count

class Peak:
    id=''
    tagNum=0
    lNum=0
    rNum=0
    lbRatio=0
    rbRatio=0
    bRatio=0
    sig=0
    expR=0
    strand='.'
    distr=[]
    def __init__(self, id, tagNum, lbRatio, rbRatio, lNum, rNum, expR, distr):
        self.id=id
        self.tagNum=tagNum
        self.lbRatio=lbRatio
        self.rbRatio=rbRatio
#         self.bRatio=abs(lbRatio-rbRatio)
        self.bRatio=lbRatio+rbRatio
        self.lNum=lNum
        self.rNum=rNum
        self.expR=expR
        self.distr=distr[:]

def zeroLog(x):
    if x<eps:
        x=eps
    return np.log(x)

def inv_logit(p):
    t1=np.exp(p)
    return t1/(1+t1)

def overlap(cent, st, end):
    return cent>st and cent<end

def readFile(infile, allPeak):
    # pr: [lowB:[tag1, tag2, tag3], medianB:[...], hiB:[...]]
    # threshold: [bRatio:[th1, th2, max], tag:[th1, th2, max]]
    nowId=''
    tagDistr=[0]*binN
    tagOne=[0]*len(tagDistr)
    tagNum=0
    lbNum=0
    rbNum=0
    lNum=0
    rNum=0
    expR=0
    atacSig=0
    fi=open(infile)
    for line in fi:
        line=line.strip()
        if len(line)<1:
            continue
        temp=line.split('\t')
        fragLen=int(temp[-2])
        b=(fragLen+offset)/binS
        if b<0 or b>=binN:
            continue
        tagDistr[b]+=1
#         if fragLen<0:
#             fragLen=-fragLen
        id='\t'.join([temp[0],temp[1],temp[2]])
        if nowId=='':
            nowId=id
            expR=float(temp[expC])
            atacSig=float(temp[atacC])
        if id != nowId:
#             print tagOne
#             if tagNum>0:
#                 for i in range(len(tagOne)):
#                     tagOne[i]/=1.0*tagNum
#                 tagOne=list(np.cumsum(tagOne))
#             allPeak.append(Peak(nowId, tagNum, 1.0*bNum/tagNum,expR))
#             allPeak.append(Peak(nowId, tagNum, 1.0*lbNum/tagNum, 1.0*rbNum/tagNum, lNum, rNum, expR))
            allPeak.append(Peak(nowId, atacSig, 1.0*lbNum/tagNum, 1.0*rbNum/tagNum, lNum, rNum, expR, tagDistr))
            nowId=id
            expR=float(temp[expC])
            atacSig=float(temp[atacC])
            tagNum=0
#             bNum=0
            lbNum=0
            rbNum=0
            lNum=0
            rNum=0
            tagDistr=[0]*binN
        if overlap((int(temp[1])+int(temp[2]))/2, int(temp[-5]), int(temp[-4])):
            if fragLen>=minB and fragLen<=maxB:
                rbNum+=1
            elif fragLen>=-maxB and fragLen<=-minB:
                lbNum+=1
        if fragLen<0:
            lNum+=1
        else:
            rNum+=1
        tagNum+=1
    fi.close()
#     allPeak.append(Peak(nowId, tagNum, 1.0*bNum/tagNum, expR))
#     allPeak.append(Peak(nowId, tagNum, 1.0*lbNum/tagNum, 1.0*rbNum/tagNum, lNum, rNum, expR))
    allPeak.append(Peak(nowId, atacSig, 1.0*lbNum/tagNum, 1.0*rbNum/tagNum, lNum, rNum, expR, tagDistr))
    if debug:
        for x in allPeak:
#             print '\t'.join([str(y) for y in [x.expR, x.tagNum, x.lbRatio, x.rbRatio, x.lbRatio+x.rbRatio]])
            print '\t'.join([str(y) for y in [x.tagNum, x.lbRatio+x.rbRatio, 1, x.expR]])

def normPk(pk, alpha):
    allTag=[x.tagNum for x in pk]
    maNum=max(allTag)
    miNum=min(allTag)
    
    allB=[x.bRatio for x in pk]
    maB=max(allB)
    miB=min(allB)
    
    for x in pk:
        if maNum-miNum<eps:
            x.tagNum=0
        else:
            x.tagNum=1.0*(x.tagNum-miNum)/(maNum-miNum)
        if maB-miB<eps:
            x.bRatio=0
        else:
            x.bRatio=1.0*(x.bRatio-miB)/(maB-miB)
        x.sig=alpha*x.tagNum+(1-alpha)*x.bRatio
#         x.sig=alpha*x.tagNum-(1-alpha)*x.bRatio
#         if x.lNum>1.5*x.rNum:
#             x.strand='+'
#         elif x.rNum>1.5*x.lNum:
#             x.strand='-'            
        if x.lbRatio>1.5*x.rbRatio:
            x.strand='+'
        elif x.rbRatio>1.5*x.lbRatio:
            x.strand='-'            

def printOut(outfile, pk):
    fo=open(outfile,'w')
    for x in pk:
        fo.write('\t'.join([str(y) for y in x.distr])+'\n')
    fo.close()

readFile(infile, allPeak)
# normPk(allPeak, alpha)
printOut(outfile, allPeak)
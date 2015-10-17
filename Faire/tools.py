#!/usr/bin/env python
from wig import Wig

from numpy import *
from wigs import Wigs
from summits import Summits

import sys,os
import math#
from copy import deepcopy
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
def aroundSummits(summits,wg,outfile,flank=1000):
    #calculate average wiggle density along the flanking regions of a set of summits
    outf=open(outfile,"w")
    step=wg.step
    flank/=step
    lst=resize([0.0],flank*2)
    num=0
    for cr in summits:
        if not wg.data.has_key(cr):continue
        lth=len(wg.data[cr])
        for smt in summits[cr]:
            smt1=smt
            smt=smt/step
            start=smt-flank
            end=smt+flank
            if start >=0:
                if end<=lth:
                    num+=1
                    for pos in range(start,end):
                        lst[pos-start]+=wg.data[cr][pos]
    for i in range(0,flank*2):
        lst[i]/=num*1.0
        outf.write(str((i-flank)*step)+"\t"+str(lst[i])+"\n")

def aroundTTS(wg,outfile,file='/pillar_storage/pillar00/kaifuc/sharon/expression/ensembl/Ensembl.tss_corrected.xls',flank=3000,title_line=True):
    #calculate average wiggle density along the flanking regions of Transcription Start Sites
    step=wg.step
    flank/=step
    #print step,flank
    outf=open(outfile,"w")
    inf=open(file)
    if title_line:title=inf.readline()
    lst=resize([0.0],flank*2)
    num=0
    for line in inf:
        col=line.split()
        chr=col[1]
        if not wg.data.has_key(chr):continue
        num+=1
        if col[2]=='+':
            tss=int(col[4])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank+i]+=wg.data[chr][tss+i]
                lst[flank-i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss-flank]
        else:
            tss=int(col[3])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank-i]+=wg.data[chr][tss+i]
                lst[flank+i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss+flank]
    lst/=num
    print num
    for i in range(0,flank*2):outf.write(str((i-flank)*step)+"\t"+str(lst[i])+"\n")
    
def aroundTSS(wg,outfile,file='/pillar_storage/pillar00/kaifuc/sharon/expression/ensembl/Ensembl.tss_corrected.xls',flank=3000,title_line=True):
    #calculate average wiggle density along the flanking regions of Transcription Terminal Sites
    step=wg.step
    flank/=step
    #print step,flank
    outf=open(outfile,"w")
    inf=open(file)
    if title_line:title=inf.readline()
    lst=resize([0.0],flank*2)
    num=0
    for line in inf:
        col=line.split()
        #col=col[1:]################### by Kaifu on June 12, 2011, to correct gene table ##########
        chr=col[1]
        if not wg.data.has_key(chr):continue
        num+=1
        if col[2]=='+':
            tss=int(col[3])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank+i]+=wg.data[chr][tss+i]
                lst[flank-i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss-flank]
        
        else:
            tss=int(col[4])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank-i]+=wg.data[chr][tss+i]
                lst[flank+i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss+flank]
        
    lst/=num # by Kaifu on June 12, 2011
    print num
    for i in range(0,flank*2):outf.write(str((i-flank)*step)+"\t"+str(lst[i])+"\n")
    

def clearBed(bed,chrsizes='/pillar_storage/pillar01/kaifuc/ref_data/sacCer1.chrom.sizes'):
    #remove the data that exceed the the ends of each chromosome
    cd={}
    for line in open(chrsizes):
        col=line.split()
        if len(col)>1:cd[col[0]]=int(col[1])
    for line in open(bed):
        col=line.split()
        if len(col)>3:
            try: int(col[2])
            except:continue
            if int(col[2])<cd[col[0]] and int(col[1])>0:print line[:-1]
def clearWig(wg):
    for line in open(wg):
        if line[0:5]=='track':continue
        else:print line[:-1]
def diffNum(file,lpv=5,pkv=50,top=100000000):
    print 'shift:',
    smt=readDiffSmt(file,category='shf',gl='',fuzzyd=40,shiftd=40,pvalue=lpv,peakv=pkv,top=top)
    print 'gain:',
    smt=readDiffSmt(file,category='dep',gl='gain',fuzzyd=40,shiftd=40,pvalue=lpv,peakv=pkv,top=top)
    print 'loss:',
    smt=readDiffSmt(file,category='dep',gl='loss',fuzzyd=40,shiftd=40,pvalue=lpv,peakv=pkv,top=top)
def fetchsmtCountArdTSS(smts,gfile,ofile=None,up=500,dn=500,step=10):
    #calculate the percentage of summits distributed in falnking region of TTS
    up-=up%step
    dn=dn-(dn%step)+step
    tsmtdic=smts.data
    smtdic={}
    for chr in tsmtdic.keys():
        smtdic[chr]={}
        for pos in tsmtdic[chr]:
            tpos=pos-pos%step
            smtdic[chr][tpos]=tsmtdic[chr][pos]
    dis={}
    num=0
    for line in open(gfile):
        if line[0]=='#':continue
        col=line.split()
        stra=col[2]
        tss=int(col[3])
        tss-=tss%step
        fr=tss-up
        to=tss+dn
        if stra!='+':
            tss=int(col[4])
            tss-=tss%step
            fr=tss-dn
            to=tss+up
        #tss-=tss%step
        #fr-=fr%step
        #to=to-(to%step)+step
        chr=col[1]
        if not dis.has_key(chr):dis[chr]={}
        if not smtdic.has_key(chr):continue
        for pos in range(fr,to,step):
            if smtdic[chr].has_key(pos):dis[chr][pos]=smtdic[chr][pos]
    num=0
    for chr in dis:num+=len(dis[chr].keys())
    print 'test',num,
    num=0
    for chr in dis:num+=len(smtdic[chr].keys())
    print num
    tre=Summits()
    tre.data=dis
    if ofile: tre.save(ofile)
    return tre

def flactdic(n):
    print "calculating flactorial for "+str(n)
    tail=1
    dic=[0]
    ftr=0
    while(tail<=n):
        ftr+=math.log(tail)
        dic.append(ftr)
        tail+=1
    return dic
def heatMapSummits(summits,wg,outfile,strd="+",flank=200,bw=0):
    outf=open(outfile,"w")
    odic={}
    step=wg.step
    flank/=step
    bw/=step
    outf.write("nuc")
    for i in range(0,flank*2):outf.write("\tpos"+str(i))
    outf.write("\n")
    #summits=summits.data
    tsmts={}
    for cr in summits:
        for smt in summits[cr]:
            tsmts[cr+" "+str(smt)]=abs(summits[cr][smt])
    from operator import itemgetter
    for item in sorted(tsmts.items(), key=itemgetter(1),reverse=True):
        #print item[1]
        #for cr in summits:
        tcol=item[0].split()
        cr,smt=tcol[0],int(tcol[1])
        if not wg.data.has_key(cr):continue
        lth=len(wg.data[cr])
        #for smt in summits[cr]:
        lst=[cr+"_"+str(smt)]
        smt1=smt
        smt=smt/step
        start=smt-flank
        end=smt+flank+1
        if start >=0:
            if end<=lth:
                for pos in range(start,end):
                    lst.append(str(wg.data[cr][pos]))
        lst=smooth(lst[1:],bw)
        if len(lst)>=flank*2 and lst[-1]!='':
            if summits[cr][smt1]<0:
                lst.reverse()
                odic[cr+"_"+str(smt1)]=lst[:-1]
            else:odic[cr+"_"+str(smt1)]=lst[1:]
            outf.write(cr+"_"+str(smt)+"\t"+"\t".join(lst)+"\n")
        else:print cr+"_"+str(smt)+"\t"+"\t".join(lst)
    return odic
def juncRatio(file,wg,ofile,flank=200):
    #juncRatio('/pillar_storage/pillar01/kaifuc/ref_data/hg18.ucsc_Genes.knownGene.table.xls',wg,'test.Qnor.xls')
    outf=open(ofile,'w')
    sum=wg.sum()
    lines=open(file).readlines()
    lines=lines[1:]
    r.options(warn=-1)
    flank/=wg.step
    for line in lines:
        col=line.split()
        if not wg.data.has_key(col[1]):continue
        juncs=col[8][:-1].split(',')
        ejuncs=col[9][:-1].split(',')
        id,eid=1,len(juncs)-1
        for i in range(id,eid):#junc in juncs:
            junc=int(juncs[i])/wg.step
            ejunc=int(ejuncs[i])/wg.step
            sjunc=int(ejuncs[i-1])/wg.step
            if junc-sjunc<flank or ejunc-junc<flank:continue
            
            if col[2]=='+':exs,exe,its,ite=junc,junc+flank,junc-flank,junc
            else:exs,exe,its,ite=junc-flank,junc,junc,junc+flank
            try:
                esum,isum=numpy.sum(wg.data[col[1]][exs:exe])/flank,numpy.sum(wg.data[col[1]][its:ite])/flank
                #ratio=esum/isum
                vec=r.c(esum,isum,sum-esum,sum-isum)#.data[chr][p],sum1-self.data[chr][p],sum2-wig2.data[chr][p])
                test=r['chisq.test'](r.matrix(vec,nrow=2))
                pvl=float(str(test[2]).split()[-1])
                if esum==0 and isum==0: pvl=1.0
                #if pvl<1e-5:
                outf.write('\t'.join([col[0],col[1],col[2],col[3],col[4],str(junc*wg.step),str(esum),str(isum),str(pvl)])+"\n")
            except:continue
        juncs=col[9][:-1].split(',')
        ejuncs=col[8][:-1].split(',')
        id,eid=1,len(juncs)-1
        for i in range(id,eid):#junc in juncs:
            junc=int(juncs[i])/wg.step
            ejunc=int(ejuncs[i])/wg.step
            sjunc=int(ejuncs[i+1])/wg.step
            if junc-ejunc<flank or sjunc-junc<flank:continue
            
            if col[2]=='+':exs,exe,its,ite=junc-flank,junc,junc,junc+flank
            else:exs,exe,its,ite=junc,junc+flank,junc-flank,junc
            esum,isum=numpy.sum(wg.data[col[1]][exs:exe])/flank,numpy.sum(wg.data[col[1]][its:ite])/flank
            #ratio=esum/isum
            try:
                vec=r.c(esum,isum,sum-esum,sum-isum)#.data[chr][p],sum1-self.data[chr][p],sum2-wig2.data[chr][p])
                test=r['chisq.test'](r.matrix(vec,nrow=2))
                pvl=float(str(test[2]).split()[-1])
                if esum==0 and isum==0: pvl=1.0
                #if pvl<1e-5:
                outf.write('\t'.join([col[0],col[1],col[2],col[3],col[4],str(junc*wg.step),str(esum),str(isum),str(pvl)])+"\n")
            except:continue
    

def lpvl(n,n1,n2,m,fd):
    minover=max(n1+n2-n,0)
    pvl=fd[n1]+fd[n-n1]+fd[n2]+fd[n-n2]-fd[n]-fd[m]-fd[n1-m]-fd[n2-m]-fd[n-n1-n2+m]
    rg=range(minover,m)
    while(len(rg)):
        over=rg.pop()
        prb=fd[n1]+fd[n-n1]+fd[n2]+fd[n-n2]-fd[n]-fd[over]-fd[n1-over]-fd[n2-over]-fd[n-n1-n2+over]
        if((pvl-prb)>700):
            return pvl #pvl=pvl
        else:
            t=pvl
            pvl=math.log(math.exp(pvl-prb)+1)+prb
            if((pvl-t)==0):return pvl
    return pvl/log(10)

def readDiffSmt(file,category='',gl='',fuzzyd=40,shiftd=40,pvalue=10,peakv=250,top=3000):
    f=open(file)
    f.readline()
    tttnum={}
    for line in f:#open(file):
        col=line[:-1].split()
        chr=col[0]
        tcol=col[5].split(":")
        t=tcol[0][:1]
        tcol[0]=int(tcol[0][1:])
        tcol[1]=int(tcol[1][1:])
        col[4]=float(col[4])
        if col[4]>=pvalue:
            if t=='c' or t=='r':
                pr=int(col[3])+tcol[0]
                pt=int(col[3])+tcol[0]+tcol[1]
            elif t=='t':
                pt=int(col[3])+tcol[0]
                pr=int(col[3])+tcol[0]+tcol[1]
            pf=int(col[3])
            pd=0
            ps=0
            if gl=='gain':
                if float(col[-2])<float(col[-1]):pd=pt
            elif gl=='left' and pr>=pt:ps=pr
            elif gl=='loss' and float(col[-2])>float(col[-1]):pd=pr
            elif gl=='right' and pr<=pt:ps=pt
            else:
                if float(col[-2])<float(col[-1]):pd=pt
                else:pd=(pr+pt)/2
                ps=pr
            if category=='shf':
                if (abs(tcol[1])>=shiftd) and abs(tcol[1])<shiftd*3:#shiftd:
                    if min(float(col[-3]),float(col[-4]))>=peakv:
                        tttnum[chr+" "+str((pr+pt)/2)]=col[4]#pt-pr
                        #print tcol[1],pt-pr
            if category=='dep':
                if True:#(abs(tcol[1])<shiftd/2):
                    if max(float(col[-1]),float(col[-2]))>=peakv:
                        if gl=='gain' and float(col[-1])>float(col[-2]):
                            if abs(float(col[4]))>=pvalue:
                                #tttnum[chr+" "+str(pf)]=col[4]#tttnum[pf]=1#
                                tttnum[chr+" "+str(pr)]=col[4]#tttnum[pf]=1#
                        if gl=='loss' and float(col[-1])<float(col[-2]):
                            if abs(float(col[4]))>=pvalue:
                                #tttnum[chr+" "+str(pf)]=col[4]#tttnum[pf]=1#
                                tttnum[chr+" "+str(pr)]=col[4]#tttnum[pf]=1#
    num=0
    print len(tttnum.keys()),
    from operator import itemgetter
    d={}
    for item in sorted(tttnum.items(), key=itemgetter(1),reverse=True):
        num+=1
        if num>top:
            smts=Summits()
            smts.data=d
            print num-1
            return smts
        col=item[0].split()
        if not d.has_key(col[0]):d[col[0]]={}
        d[col[0]][int(col[1])]=float(item[1])
    smts=Summits()
    smts.data=d
    print num-1
    return smts


def rpvl(n,n1,n2,m,fd):
    maxover=min(n1,n2)
    pvl=fd[n1]+fd[n-n1]+fd[n2]+fd[n-n2]-fd[n]-fd[m]-fd[n1-m]-fd[n2-m]-fd[n-n1-n2+m]
    for over in range((m+1),(maxover+1)):
        prb=fd[n1]+fd[n-n1]+fd[n2]+fd[n-n2]-fd[n]-fd[over]-fd[n1-over]-fd[n2-over]-fd[n-n1-n2+over]
        if((pvl-prb)>700):
            return pvl #pvl=pvl
        else:
            t=pvl
            pvl=math.log(math.exp(pvl-prb)+1)+prb
            if((pvl-t)==0):
                return pvl
    return pvl/log(10)

def readPeakSmt(file,minv=1,maxv=500000000,top=0,bottom=100000):
    print 'reading summits from',file
    f=open(file)
    f.readline()
    summit={}
    count={}
    n=0
    tttnum={}
    for line in f:#open(file):
        col=line[:-1].split()
        if not summit.has_key(col[0]):summit[col[0]]={}
        #ttt=col[0].split("_")
        #if not summit.has_key(ttt[0]):summit[ttt[0]]={}
        if float(col[4])<=maxv:
            if float(col[4])>=minv:
                n+=1
                #summit[col[0]][(int(col[1])+int(col[2]))/2]=float(col[4])
                tttnum[col[0]+" "+col[3]]=float(col[4])
                #summit[ttt[0]][int(ttt[1])]=float(col[4])
    
    num=0
    from operator import itemgetter
    d={}
    #stop=top
    items=sorted(tttnum.items(), key=itemgetter(1),reverse=True)
    #if bottom>0:
    #    items=sorted(tttnum.items(), key=itemgetter(1))
    #    stop=bottom
    for item in items:
        num+=1
        if num>top:
            col=item[0].split()
            if not d.has_key(col[0]):d[col[0]]={}
            d[col[0]][int(col[1])]=float(item[1])
        if num>=bottom:
            #smts=Summits()
            #smts.data=d
            print 'top',top+1,'to',num,'of',len(tttnum.keys())
            print 'cutoff:',item[1]
            return d
    #smts=Summits()
    #smts.data=d
    print 'top',top+1,'to',num,'of',len(tttnum.keys())
    return d

def readStableSmt(file,minv=0,maxv=500000000):
    print 'reading summits from',file
    f=open(file)
    f.readline()
    summit={}
    count={}
    n=0
    for line in f:#open(file):
        col=line[:-1].split()
        if not summit.has_key(col[0]):summit[col[0]]={}
        #ttt=col[0].split("_")
        #if not summit.has_key(ttt[0]):summit[ttt[0]]={}
        col[-1],col[-2]=float(col[-1]),float(col[-2])
        fd=min(col[-1],col[-2])/max(col[-1],col[-2])
        if fd>1/1.1:
            if max(float(col[-1]),float(col[-2]))<=maxv:
                if min(float(col[-1]),float(col[-2]))>=minv:
                    n+=1
                    summit[col[0]][(int(col[1])+int(col[2]))/2]=max(float(col[-1]),float(col[-2]))
                    #summit[ttt[0]][int(ttt[1])]=float(col[4])
    print n
    smt=Summits()
    smt.data=summit
    return smt
    


def scatterSummits(summits,wggc,wg0,wg1,wg2,outfile,flank=40,bin=5):
    outf=open(outfile+'.xls',"w")
    step=wg0.step
    if step!=wg1.step:return 0
    if step!=wg2.step:return 0
    flank/=step
    sdic={}
    for chr in summits:
        lth=numpy.size(wg0.data[chr])
        for smt in summits[chr]:
            smt/=step
            for pos in range(smt-flank,smt+flank+1):
                if pos>0 and pos<lth:sdic[str(chr)+"\t"+str(pos)]=wg2.data[chr][pos]#numpy.sum(wg0.data[chr][(pos-5):(pos+6)])#wg0.data[chr][pos]
    slst=sorted([(value,key) for (key,value) in sdic.items()])
    for i in range(0,len(slst)-bin,bin):
        vgc,v0,vd,v1,v2=0,0,0,0,0
        for (value,key) in slst[i:(i+bin)]:
            chr,smt=key.split()
            smt=int(smt)
            vgc+=numpy.sum(wggc.data[chr][(smt-5):(smt+6)])
            v0+=wg0.data[chr][smt]
            #v0+=numpy.sum(wg0.data[chr][(smt-5):(smt+6)])
            vd+=wg1.data[chr][smt]-wg2.data[chr][smt]
            v1+=wg1.data[chr][smt]
            v2+=wg2.data[chr][smt]
        vgc=vgc*10/11
        outf.write(str(i)+"\t"+str(vgc/bin)+"\t"+str(v0/bin)+"\t"+str(vd/bin)+"\t"+str(v1/bin)+"\t"+str(v2/bin)+"\n")


def smooth(lst,bw=0):
    if bw==0:return lst
    lth=len(lst)
    for i in range(0,lth):lst[i]=float(lst[i])
    t=lst[(0-bw):]+lst+lst[:bw]
    for i in range(bw,bw+lth):
        for j in range(1,bw+1):lst[i-bw]+=t[i-j]+t[i+j]
        lst[i-bw]/=bw*2+1
    for i in range(0,lth):lst[i]=str(lst[i])
    return lst    
def stablePeak(rsmt,tsmt,dwig,rwig,twig,file,step=10,dis=150):
    dic={}
    outf=open(file,"w")
    #outf.write("chr\tviewStart\tviewEnd\tsummit\t-log10(P value)\tfuzzy:shiftDistance\tcontrolSummit\ttreatSummit\tcontrolPoint\ttreatPoint\n")
    for cr in rsmt.data:
        ps=rsmt.data[cr].keys()
        ps.sort()
        for p in ps:
            d1=localNeighbour(tsmt.data[cr],p,step,dis)
            #d2=0-localNeighbour(rsmt.data[cr],p+d1,step,dis)
            if d1==0:
                pv=0
                p=p/step
                if p>0 and p<numpy.size(dwig.data[cr])and p<numpy.size(rwig.data[cr])and p<numpy.size(twig.data[cr]):pv=abs(dwig.data[cr][p])
                if pv<0.5:
                    if twig.data[cr][p]>50 and twig.data[cr][p]>50:
                        if not dic.has_key(cr):dic[cr]={}
                        if not dic[cr].has_key(p):dic[cr][p]=1
                        outf.write("\t".join([cr,str(p-500),str(p+500),str(p*step),str(0),'+',str(pv),str(d1),str(twig.data[cr][p]),str(rwig.data[cr][p])])+"\n")
    smt=Summits()
    smt.data=dic
    return smt

def summitDis(summits,bgsmts,markfile,ofile,gsize=12070970,step=10):
    #the statistics distribution of summits in different genomic regions relative to the gene structure.
    outf=open(ofile,'w')
    gsize=gsize/step
    chr=""
    dic={}
    bgdic={}
    for k in ['inter','TSSup','TSSdn','exon','TTSup','TTSdn','all']:
        dic[k]=0
        bgdic[k]=0
    #dic={'all':0,'inter':0}
    #bgdic={'all':0,'inter':0}
    #cat={'all':0,'inter':0}
    for chr in summits:
        dic['all']+=len(summits[chr].keys())
    for chr in bgsmts:
        bgdic['all']+=len(bgsmts[chr].keys())
    for line in open(markfile):
        if line[0]=='v':
            col=line.split()
            tcol=col[1].split("=")
            chr=tcol[1]
        else:
            if not summits.has_key(chr):continue
            if not bgsmts.has_key(chr):continue
            #cat['all']+=1
            col=line.split()
            pos=int(col[0])
            pos-=pos%step
            #for r in col[1:]:
            #    if cat.has_key(r):cat[r]+=1
            #    else:cat[r]=1
            if summits[chr].has_key(pos):
                for r in col[1:]:
                    if dic.has_key(r):dic[r]+=1
                    else:dic[r]=1
                summits[chr].pop(pos)
            if bgsmts[chr].has_key(pos):
                for r in col[1:]:
                    if bgdic.has_key(r):bgdic[r]+=1
                    else:bgdic[r]=1
                bgsmts[chr].pop(pos)
    for chr in summits:
        dic['inter']+=len(summits[chr].keys())
    for chr in bgsmts:
        bgdic['inter']+=len(bgsmts[chr].keys())
    #cat['inter']=gsize-cat['all']
    #cat['all']=gsize
    from math import log
    #print dic['all'],bgdic['all']
    outf.write('\tlog2FC\tlog10PV\tsmt_num\tbg_num\n')
    for k in ['TTSdn','inter','TSSup','TSSdn','exon','TTSup','all']:
        from rpy2.robjects import r
        r.options(warn=-1)
        from math import log10
        vec=r.c(bgdic[k],dic[k],bgdic['all'],dic['all'])#,sum1-self.data[chr][p],sum2-wig2.data[chr][p])
        test=r['chisq.test'](r.matrix(vec,nrow=2))
        pvl=float(str(test[2]).split()[-1])
        #if bgdic[k]*1.0/bgdic['all']>dic[k]*1.0/dic['all']:pvl=str(log10(pvl))
        #else:
        #pvl=str(-log10(pvl))
        #fc=(dic[k]*1.0/dic['all'])/(bgdic[k]*1.0/bgdic['all'])
        #outf.write(str(k)+"\t"+str(log(fc)/log(2))+"\t"+str(pvl)+"\t"+str(dic[k])+"\t"+str(bgdic[k])+"\n")
        outf.write(str(k)+"\t"+str(dic[k])+"\t"+str(bgdic[k])+"\n")
def smtdis(smtfile,gmarkfile):
    smt={}
    all=0
    for line in open(smtfile).readlines()[1:]:
        col=line.split()
        if not smt.has_key(col[0]):smt[col[0]]={}
        smt[col[0]][int(col[3])-int(col[3])%10]=1
        all+=1
    ot={}
    chr=''
    for line in open(gmarkfile):
        if line[0]=='v':
            col=line.split()
            tcol=col[1].split("=")
            chr=tcol[1]
        else:
            col=line.split()
            p=int(col[0])-1
            if smt.has_key(chr):
                if smt[chr].has_key(p):
                    if 'TSSup' in col[1:]:
                        if not ot.has_key('TSSup'):ot['TSSup']=1
                        else:ot['TSSup']+=1
                    elif 'TSSdn' in col[1:]:
                        if not ot.has_key('TSSdn'):ot['TSSdn']=1
                        else:ot['TSSdn']+=1
                    else:
                        for k in col[1:]:
                            if not ot.has_key(k):ot[k]=1
                            else:ot[k]+=1
    print 'all',all
    for k in ot:
        print k,ot[k]
        
def smtCountArdTSS(smts,gfile,ofile,flank=2500,step=10):
    #calculate the percentage of summits distributed in falnking region of TTS
    tsmtdic=smts#.data
    smtdic={}
    for chr in tsmtdic.keys():
        smtdic[chr]={}
        for pos in tsmtdic[chr]:
            tpos=pos-pos%step
            smtdic[chr][tpos]=tsmtdic[chr][pos]
    dis={}
    num=0
    for line in open(gfile):
        if line[0]=='#':continue
        num+=1
        col=line.split()
        stra=col[2]
        try: tss=int(col[3])
        except:continue
        if stra!='+':tss=int(col[4])
        tss-=tss%step
        chr=col[1]
        if not smtdic.has_key(chr):continue
        for pos in range(tss-flank,tss+flank,step):
            if smtdic[chr].has_key(pos):
                if stra=='+':
                    if dis.has_key(pos-tss):dis[pos-tss]+=1#smtdic[chr][pos]
                    else: dis[pos-tss]=1#smtdic[chr][pos]
                else:
                    if dis.has_key(tss-pos):dis[tss-pos]+=1#smtdic[chr][pos]
                    else:dis[tss-pos]=1#smtdic[chr][pos]
    num*=step
    outf=open(ofile,"w")
    for k in range(0-flank,flank+1,step):
        if not dis.has_key(k):dis[k]=0
        outf.write(str(k)+"\t"+str(dis[k]*1.0/num)+"\n")

def smtCountArdTTS(smts,gfile,ofile,flank=2500,step=10):
    #calculate the percentage of summits distributed in falnking region of TTS
    tsmtdic=smts#.data
    smtdic={}
    for chr in tsmtdic.keys():
        smtdic[chr]={}
        for pos in tsmtdic[chr]:
            tpos=pos-pos%step
            smtdic[chr][tpos]=tsmtdic[chr][pos]
    dis={}
    num=0
    for line in open(gfile):
        num+=1
        if line[0]=='#':continue
        col=line.split()
        stra=col[2]
        try:tss=int(col[4])
        except:continue
        if stra!='+':tss=int(col[3])
        tss-=tss%step
        #print tss
        chr=col[1]
        if not smtdic.has_key(chr):continue
        for pos in range(tss-flank,tss+flank,step):
            if smtdic[chr].has_key(pos):
                
                if stra=='+':
                    if dis.has_key(pos-tss):dis[pos-tss]+=1#smtdic[chr][pos]
                    else: dis[pos-tss]=1#smtdic[chr][pos]
                else:
                    if dis.has_key(tss-pos):dis[tss-pos]+=1#smtdic[chr][pos]
                    else:dis[tss-pos]=1#smtdic[chr][pos]
    num*=10
    #for chr in smtdic:num+=len(smtdic[chr].keys())
    outf=open(ofile,"w")
    for k in range(0-flank,flank+1,step):
        if not dis.has_key(k):dis[k]=0
        outf.write(str(k)+"\t"+str(dis[k]*1.0/num)+"\n")

def smtIntesect(smts1,smts2,dis=0,file=None):
    s1=smts1.data#readPeakSmt(f1).data
    s2=smts2.data#readPeakSmt(f2).data
    n,n1,n2=0,0,0
    com={}
    spe1={}
    spe2={}
    for k in s1:
        if not com.has_key(k):com[k]={}
        if not spe1.has_key(k):spe1[k]={}
        if not spe2.has_key(k):spe2[k]={}
        for p in s1[k].keys():
            n1+=1
            a=0
            for i in range(p-dis,p+dis+1):
                if s2[k].has_key(i):a+=1
            if a:
                com[k][p]=1
                n+=1
            else:spe1[k][p]=1
    for k in s2:
        if not com.has_key(k):com[k]={}
        if not spe1.has_key(k):spe1[k]={}
        if not spe2.has_key(k):spe2[k]={}
        for p in s2[k].keys():
            n2+=1
            for i in range(p-dis,p+dis+1):
                if s1[k].has_key(i):a+=1
            if a:
                com[k][p]=1
            else:spe2[k][p]=1
    out=Summits()
    out.data=spe1
    if file:out.save(file,append=True)
    out.data=spe2
    if file:out.save(file,append=False)
    out.data=com
    if file:out.save(file,append=True)
    print 'group1',n1,n1-n
    print 'group2',n2,n2-n
    print 'union:',n1+n2-n
    print 'overlap:',n
    print 'Jaccard index:',n*1.0/(n1+n2-n)
    fd=flactdic(64456)
    try:
        print rpvl(64456,n1,n2,n,fd),10**rpvl(64456,n1,n2,n,fd),"\n\n"
    except:
        return 0    

    return out
def smtVenndiagram(f1,f2,f3,dis=40):
    s1=readPeakSmt(f1).data
    s2=readPeakSmt(f2).data
    s3=readPeakSmt(f3).data
    n12,n13,n23,n123=0,0,0,0
    n1,n2,n3=0,0,0
    for k in s1:
        for p in s1[k].keys():
            n1+=1
            a12,a13=0,0
            for j in range(0,dis+1):
                #for i in range(p-dis,p+dis+1):
                if a12==0:
                    if s2[k].has_key(p-j):a12=p-j
                    if s2[k].has_key(p+j):a12=p+j
                if a13==0:
                    if s3[k].has_key(p-j):a13=p-j
                    if s3[k].has_key(p+j):a13=p+j
            if a12: n12+=1
            if a13: n13+=1
            if abs(a12-a13)<=dis:n123+=1
    for k in s2:
        for p in s2[k].keys():
            n2+=1
            a23=0
            for j in range(0,dis+1):
                    if s3[k].has_key(p-j):a23=p-j
                    if s3[k].has_key(p+j):a23=p+j
            if a23: n23+=1
    for k in s3:
        for p in s3[k].keys():n3+=1
    print 'intersections:'
    print 'all',n123
    print f1,f2,n12-n123
    print f1,f3,n13-n123
    print f2,f3,n23-n123
    print f1,n1-n12-n13+n123
    print f2,n2-n23-n12+n123
    print f3,n3-n13-n23+n123
    print 'union:',n123+n12-n123+n13-n123+n23-n123+n1-n12-n13+n123+n2-n23-n12+n123+n3-n13-n23+n123
    print 'Jaccard index:',n123*1.0/(n123+n12-n123+n13-n123+n23-n123+n1-n12-n13+n123+n2-n23-n12+n123+n3-n13-n23+n123)
    
def smt2tss(smt,gene_table_file='../Ensembl.xls',ofile='result.xls',up=0,dn=500,step=1,tss=0):
    #match summits to transcription start sites
    fo=open(ofile,'w')
    #smt=smt.data
    gf=open(gene_table_file)#/pillar_storage/pillar01/kaifuc/ref_data/sacCer1_Ensemblgenes.xls 
    fo.write(gf.readline()[:-1]+"\t"+'neighbour summits\n')
    for line in gf:
        col=line.split()
        if not smt.has_key(col[1]):continue
        tup,tdn,pos=up,dn,int(col[3])+tss
        if col[2]!='+':tup,tdn,pos=dn,up,int(col[4])-tss
        ps=[]
        for p in range(pos-tup,pos+tdn,step):
            if smt[col[1]].has_key(p):ps.append(p)
        if len(ps)>0:
            for i in range(0,len(ps)):ps[i]=str(ps[i])
            fo.write(line[:-1]+"\t"+"\t".join(ps)+"\n")

def overlap3(par=[]):
    dic0={}
    for line in open(par[0]):dic0[line[:-1]]=1
    n=len(dic0.keys())
    dic1={}
    f2=open(par[1])
    for line in f2:
        col=line.split()
        if dic0.has_key(col[0]):
            dic1[col[0]]=1
    n1=len(dic1.keys())
    f3=open(par[2])
    dic2={}
    dic3={}
    for line in f3:
        col=line.split()
        if dic0.has_key(col[0]):
            dic2[col[0]]=1
            if dic1.has_key(col[0]):
                dic3[col[0]]=1
    n2=len(dic2.keys())
    n3=len(dic3.keys())
    print n1,n2,n3

    fd=flactdic(n)
    try:
        print rpvl(n,n1,n2,n3,fd),10**rpvl(n,n1,n2,n3,fd),"\n\n"
    except:
        return 0
def test1():
    from math import log10
    dic={}
    for line in open(sys.argv[1]):
        if line[:2]=='ID':continue
        col=line.split()
        if float(col[1])>0:dic[col[0]]=-log10(float(col[-2]))
        else:dic[col[0]]=-log10(float(col[-2]))
    for line in open(sys.argv[2]):
        col=line.split()
        if dic.has_key(col[0]):print col[0]+"\t"+str(dic[col[0]])+"\t"+"\t".join(col[1:])



def overlap2(par=[]):
    dic1={}
    f1=open(par[0])
    for line in f1:
        col=line.split()
        dic1[col[0]]=1
    n1=len(dic1.keys())
    f2=open(par[1])
    dic2={}
    for line in f2:
        col=line.split()
        dic2[col[0]]=1
    n2=len(dic2.keys())
    n3=0
    for k in dic2.keys():
        if dic1.has_key(k):
            n3+=1
            print k
    #fd=flactdic(17533)
    #print n1,n2,n3
    #print rpvl(17533,n1,n2,n3,fd),10**rpvl(17533,n1,n2,n3,fd),"\n\n"


def temp_aroundTSS(wg,outfile,file='/pillar_storage/pillar00/kaifuc/sharon/expression/ensembl/Ensembl.tss_corrected.xls',flank=3000,title_line=True):
    #calculate average wiggle density along the flanking regions of Transcription Terminal Sites
    step=wg.step
    flank/=step
    #print step,flank
    outf=open(outfile,"w")
    inf=open(file)
    if title_line:title=inf.readline()
    lst=resize([0.0],flank*2)
    num=0
    for line in inf:
        col=line.split()
        #col=col[1:]################### by Kaifu on June 12, 2011, to correct gene table ##########
        chr=col[0]
        if not wg.data.has_key(chr):continue
        num+=1
        if col[2]=='+':
            tss=int(col[1])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank+i]+=wg.data[chr][tss+i]
                lst[flank-i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss-flank]
        else:
            tss=int(col[1])/step
            if (tss+flank)>=wg.chrSize(chr)/step:continue
            if tss<flank:continue
            lst[flank]+=wg.data[chr][tss]
            for i in range(1,flank):
                lst[flank-i]+=wg.data[chr][tss+i]
                lst[flank+i]+=wg.data[chr][tss-i]
            lst[0]+=wg.data[chr][tss+flank]
    lst/=num # by Kaifu on June 12, 2011
    print num
    for i in range(0,flank*2):outf.write(str((i-flank)*step)+"\t"+str(lst[i])+"\n")
    
if __name__=='__main__':
    
    smtdis(smtfile='simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.xls',gmarkfile='/lilab/kaifuc/ref_data/sacCer1_Ensemblgenes.genomeMark-TSSup350TSSdn50TTSup0TTSdn0.xls')
    #smtdis(smtfile='simuDefaultDiffNuc.eth-ypd.1e-07.shft50-90.xls',gmarkfile='/lilab/kaifuc/ref_data/sacCer1_Ensemblgenes.genomeMark-TSSup350TSSdn50TTSup0TTSdn0.xls')
    #smtdis(smtfile='simuDefaultDiffNuc.eth-ypd.1e-07.fuzz1e-07.xls',gmarkfile='/lilab/kaifuc/ref_data/sacCer1_Ensemblgenes.genomeMark-TSSup350TSSdn50TTSup0TTSdn0.xls')
    '''
    oc={}
    for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.xls'):
        col=line.split()
        oc[col[0]+'_'+col[3]]=1
    sh={}
    for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.shft50-90.xls'):
        col=line.split()
        sh[col[0]+'_'+col[3]]=1
    fz={}
    for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.fuzz1e-07.xls'):
        col=line.split()
        fz[col[0]+'_'+col[3]]=1
    n12,n13,n23,n123,n1,n2,n3=0,0,0,0,0,0,0
    for k in sh:
        n1+=1
        if k in oc:n12+=1
        if k in fz:
            n13+=1
            if k in oc:n123+=1
    print n12,n13,n123
    print len(sh.keys()),len(oc),len(fz)
        
    
    smt={}
    for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part1.xls'):
        #for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.xls'):
        col=line.split()
        if not smt.has_key(col[0]):smt[col[0]]={}
        smt[col[0]][int(col[3])]=1
    #smt2tss(smt,gene_table_file='/lilab/kaifuc/ref_data/sacCer1_Ensemblgenes.xls',ofile='result.xls',up=350,dn=50,step=1,tss=0)
    wg=Wig('../../default/Eth-ypd/Eth-ypd.smooth.control.wig')
    heatMapSummits(smt,wg,outfile='heatmap_control_simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part1.xls',flank=200,bw=0)
    wg=Wig('../../default/Eth-ypd/Eth-ypd.smooth.treat.wig')
    heatMapSummits(smt,wg,outfile='heatmap_treat_simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part1.xls',flank=200,bw=0)
    
    smt={}
    for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part2.xls'):
        #for line in open('simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.xls'):
        col=line.split()
        if not smt.has_key(col[0]):smt[col[0]]={}
        smt[col[0]][int(col[3])]=1
    #smt2tss(smt,gene_table_file='/lilab/kaifuc/ref_data/sacCer1_Ensemblgenes.xls',ofile='result.xls',up=350,dn=50,step=1,tss=0)
    wg=Wig('../../default/Eth-ypd/Eth-ypd.smooth.control.wig')
    heatMapSummits(smt,wg,outfile='heatmap_control_simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part2.xls',flank=200,bw=0)
    wg=Wig('../../default/Eth-ypd/Eth-ypd.smooth.treat.wig')
    heatMapSummits(smt,wg,outfile='heatmap_treat_simuDefaultDiffNuc.eth-ypd.1e-07.occ1e-07.fuzz1e-5higher.part2.xls',flank=200,bw=0)
    
    d={}
    for line in open('simuDefaultDiffNucGenes.eth-ypd.1e-07.shft50-90.xls'):d[line.split()[0]]=1
    for line in open('simuDefaultDiffNucGenes.eth-ypd.1e-07.fuzz1e-07.xls'):d[line.split()[0]]=1
    for line in open('result.xls'):
        col=line.split()
        if not d.has_key(col[0]):print col[0]
    '''
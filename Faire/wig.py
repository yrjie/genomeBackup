#!/usr/bin/env python
import numpy,re
from copy import deepcopy
from rpy2.robjects import r,FloatVector
from math import log10,sqrt,log
from summits import Summits
from time import time

class Wig:
    def __init__(self,file="",gfile='',step=0,suppress=False):
        '''
        Parameter:
            file: a path to a file in Wiggle format
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region.
        '''
        self.data = {} # a dictionary in the format self.data[chr_name]=numpy_array
        self.step=step
        if gfile!='':
            for line in open(gfile):
                col=line.split()
                self.data[col[0]]=numpy.array([0.0])
                if step>0:self.data[col[0]].resize(int(col[1])/step+1,refcheck=0)
                else:self.data[col[0]].resize(int(col[1])+1,refcheck=0)
        if file !="": self.load(file,gfile=gfile,suppress=suppress)
    def absSum(self):
        '''
        Description:
            the sum of absolute value at each data point.
        
        Parameter:
            None
            
        Value:
            float value
        '''
        wig=self.data
        sum=0.0
        for chr in wig:
            for v in wig[chr]:sum+=abs(v)
        return sum
    def add(self,wig2):
        '''
        Description:
            add value at each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance whose value will not change but will be add to the Wig class instance calling this method.
        
        Value:
            None
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        chrs={}
        for chr in self.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            if not self.data.has_key(chr):self.addChr(chr)
            if not wig2.data.has_key(chr):wig2.addChr(chr)
            lth1=self.data[chr].size
            lth2=wig2.data[chr].size
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            elif lth1<lth2:self.data[chr].resize(lth2,refcheck=0)
            self.data[chr]+=wig2.data[chr]
        return True
    def addChr(self,chr):
        '''
        Description:
            add a new chrosome to the existance Wig class instance
        
        Parameter:
            chr: the name of the chrosome to be added
            
        Value: none
        '''
        self.data[chr]=numpy.array([0.0])
    def ajust_size(self,gfile):
        step=self.step
        crs={}
        for line in open(gfile):
            col=line.split()
            size=int(col[1])/step+1
            crs[col[0]]=size
            if size<self.data[col[0]].size:self.data[col[0]]=self.data[col[0]][:size] #resize(int(col[1])/step+1,refcheck=0)
            if size>self.data[col[0]].size:self.data[col[0]].resize(size,refcheck=0)
        tcrs=self.data.keys()
        for cr in tcrs:
            if not crs.has_key(cr):self.data.pop(cr)
        
    def getChr(self,chr):
        '''
        Description:
            retrieve a chrosome by name.
            
        Parameter:
            chr: the name of the chrosome to be retrived.
        
        Value:
            a list as an instance of the numpy.array.
        '''
        return self.data[chr]
    def correlation(self,wg2,ofile,mind=0,maxd=50):
        '''
        Description:
            Calculate correlation between the two Wig class instances allowing shift distances.
            
        Parameters:
            ofile: the output file used to save the result.
            mind: minimal shift distance
            maxd: maximal shift distance
        
        Value:
            None
        '''
        d={}
        chrs=self.getChrs()
        chrs2=wg2.getChrs()
        mind/=self.step
        maxd/=self.step
        for chr in chrs:
            d[chr]={}
            if chr in chrs2:
                cs=max(self.chrSize(chr),wg2.chrSize(chr))
                if self.chrSize(chr)!=cs:self.resizeChr(chr,cs)
                if wg2.chrSize(chr)!=cs:wg2.resizeChr(chr,cs)
                for td in range(mind,maxd+1):
                    v=r.cor(FloatVector(self.getChr(chr)[:(cs/self.step-td)]),FloatVector(wg2.getChr(chr)[td:cs/self.step]))
                    d[chr][td]=float(str(v).split()[-1])
        fo=open(ofile,'w')
        fo.write('Shift_disance\tcorrelation_coefficient\n')
        for td in range(mind,maxd+1):
            v=0.0
            s=0
            gs=self.gsize()
            for chr in d:
                v+=self.chrSize(chr)*d[chr][td]
            fo.write(str(td*self.step)+'\t'+str(v/gs)+'\n')

    def bgsub(self,wig2,lmd=1000):
        '''
        Description:
            subtract the value of each data point in wig2 from self after smoothing of wig2
        
        Parameter:
            wig2: the Wig class instance used to do the subtraction
            lmd: the bin size used to smooth wig2
        
        Value:
            None
        
        Note:
            a copy of wig2 is smoothed and used to do the subtraction, wig2 will not change.
        '''
        wig2=deepcopy(wig2)
        if self.step!=wig2.step:wig2.changeStep(self.step)
        chrs={}
        for chr in self.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            lth1=len(self.data[chr])
            lth2=len(wig2.data[chr])
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            else:self.data[chr].resize(lth2,refcheck=0)
        if lmd>0:wig2.smooth(lmd=lmd)
        wig2.foldChange(self.sum()/wig2.sum())
        print 'before subtracting:',self.sum()
        self.subtract(wig2)
        print 'after subtracting:',self.sum()
        self.rvNeg()
        print 'after removing negtive values:',self.sum()
        return True
    def callRegions(self,ofile=None,width=0,distance=165,pheight=0,height=0,calculate_P_value=1,mode='w',title_line=1,pos_only=False,fold=0,suppress=False,fdr=False,fdrSampleSize=1000000):
        '''
        Description:
            This fuction is designed to call broad peaks, such as histone modification peaks.
        
        Parameter:
            ofile: a path to the file used to save the peaks.
            width: minimal width of peaks
            distance: minimal distance between peaks, neighboring peaks with distance shorter than this value will be merged.
            pheight: a P value cutoff used to call peaks.
            height: the occupancy cutoff used to call peaks. valide only when pheight is set to 1.
            calculate_P_value: calculate P value for each peak if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
        
        Value:
            pks[chrosome_name][start_position]=end_position
            
        '''
        ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE))}''')
        m=self.mean()
        if height==0 and pheight!=0:
            #if fold!=0:height=m*fold
            #else:
            height=int(m+0.5)
            while( (0-float(str(ppois(height,m)).split()[-1])/log(10)) < pheight):height+=1
        #if not suppress:
        print 'whole genome aveage value is '+str(m)+', use calling cutoff',height
        dic=self.data
        step=self.step
        twidth=width/step
        if twidth<1:twidth=1
        pks={}
        if not suppress: print 'calling ...'
        for chr in dic:
            pks[chr]={}
            lth=dic[chr].size#.chrSize(chr)/twstep-twidth
            if lth<=0:continue
            if not suppress: print chr,
            start,pos=-1,0
            while pos <lth:
                if start<0:
                    if dic[chr][pos]>=height:start=pos
                elif dic[chr][pos]<height:
                    if pos-start>=twidth:pks[chr][start*step]=(pos)*step
                    start=-1
                pos+=1
            if not suppress: print len(pks[chr])
        if distance>0: 
            if not suppress: print 'mering'
            crs=pks.keys()
            crs.sort()
            for chr in crs :
                ps=pks[chr].keys()
                if not suppress: print chr, 'from',len(ps),
                ps.sort()
                lth=len(ps)-1
                i=0
                while i<lth:
                    #j=i+1
                    #while ps[j]-pks[chr][ps[i]]<=distance:
                    if ps[i+1]-pks[chr][ps[i]]<=distance:# or pks[chr][ps[i+1]]<=pks[chr][ps[i]]:
                        if pks[chr][ps[i]]<pks[chr][ps[i+1]]:pks[chr][ps[i]]=pks[chr][ps[i+1]]
                        pks[chr].pop(ps[i+1])
                        ps[i+1]=ps[i]
                    i+=1
                if not suppress: print 'to',len(pks[chr])
        if ofile!=None:
            outf=open(ofile,mode)
            if title_line:
                if mode=='w':
                    if calculate_P_value:outf.write("chr\tstart\tend\tsummit_pos\tsummit_value\tstrand\ttotal_signal\twidth_above_cutoff\tsummit_minus_logP\n")
                    else:outf.write("chr\tstart\tend\tsummit_pos\tsummit_value\tstrand\ttotal_signal\twidth_above_cutoff\n")
            crs=pks.keys()
            crs.sort()
            for chr in crs :
                #out[chr]={}
                ps=pks[chr].keys()
                ps.sort()
                for p in ps:
                    width_above_cutoff=0
                    s,e=p/step,pks[chr][p]/step
                    v=max(dic[chr][s:e])
                    auc=0#dic[chr][s:e].sum()*step
                    if ((e-s)>=twidth):# and (v>=sheight):
                        smts=[]
                        for i in range(s,e):
                            if dic[chr][i]==v:smts.append(str(i*step))
                            if dic[chr][i]>=height:width_above_cutoff+=step
                            if pos_only:
                                if dic[chr][i]>0:auc+=dic[chr][i]
                            else:auc+=dic[chr][i]
                        auc=auc*step
                        smt=','.join(smts)
                        lth=len(smts)/2
                        if calculate_P_value:
                            pvl=float(str(ppois(v,m)).split()[-1])/log(10)
                            outf.write(chr+"\t"+str(p)+"\t"+str(pks[chr][p])+"\t"+smt+'\t'+str(v)+'\t+\t'+str(auc)+'\t'+str(width_above_cutoff)+'\t'+str(0-pvl)+"\n")
                        else:outf.write(chr+"\t"+str(p)+"\t"+str(pks[chr][p])+"\t"+smt+'\t'+str(v)+'\t+\t'+str(auc)+'\t'+str(width_above_cutoff)+'\n')
                    else:pks[chr].pop(p)
        return pks
        
    def fillRegions(self,regions={},file=None,pheight=1,height=0,width=0,calculate_P_value=1,pos_only=False,suppress=False):
        '''
        Add description here
        '''
        ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE))}''')
        m=self.mean()
        if height==0 and pheight!=0:
            #if fold!=0:height=m*fold
            #else:
            height=int(m+0.5)
            while( (0-float(str(ppois(height,m)).split()[-1])/log(10)) < pheight):height+=1
        #if not suppress:
        print 'whole genome aveage value is '+str(m)+', use cutoff',height
        #lines=open(file).readlines()
        outf=open(file,'w')
        if calculate_P_value==1:outf.write('chr\tstart\tend\tcenter\twidth_above_cutoff\ttotal_signal\theight\theight_logP\n')
        else:outf.write('chr\tstart\tend\tcenter\twidth_above_cutoff\ttotal_signal\theight\n')
        step=self.step
        dic=self.data
        total_width_above_cutoff=0
        #for line in lines[1:]:
        #    col=line.split()
        for chr in regions:
            starts=regions[chr].keys()
            starts.sort()
            for start in starts:
                if regions[chr][start]-start<width:continue
                s,e=start/step,regions[chr][start]/step+1
                if dic[chr].size<=e:dic[chr].resize(e+1,refcheck=0)
                width_above_cutoff=0
                v=max(dic[chr][s:e])
                if pos_only:v=max(v,0)
                auc=0
                #substart=-1
                #subpeaks=''
                i=s
                while i<e:
                    '''
                    if substart<0:
                        if dic[chr][i]>=height:substart=i
                    elif dic[chr][i]<height:
                        subpeaks+=str(substart*step)+':'+str(i*step)+','
                        substart=-1
                    '''
                    if dic[chr][i]>=height:width_above_cutoff+=step
                    if pos_only:
                        if dic[chr][i]>0:auc+=dic[chr][i]
                    else:auc+=dic[chr][i]
                    i+=1
                #if substart>=0:subpeaks+=str(substart*step)+':'+str((e-1)*step)+','
                auc=auc*step
                total_width_above_cutoff+=width_above_cutoff
                if calculate_P_value:
                    pvl=float(str(ppois(v,m)).split()[-1])/log(10)
                    outf.write(chr+'\t'+str(start)+'\t'+str(regions[chr][start])+'\t'+str((regions[chr][start]+start)/2)+'\t'+str(width_above_cutoff)+'\t'+str(auc)+'\t'+str(v)+'\t'+str(0-pvl)+"\n")
                else:outf.write(chr+'\t'+str(start)+'\t'+str(regions[chr][start])+'\t'+str((regions[chr][start]+start)/2)+'\t'+str(width_above_cutoff)+'\t'+str(auc)+'\t'+str(v)+'\n')
        print 'total_width_above_cutoff:',total_width_above_cutoff
        outf.close()
        
    def callPositions(self,ofile,width=40,distance=165,edge=1,pcut=1e-5,height=0,fill_gap=False,fill_value=1,calculate_P_value=1,mode='w',title_line=1,poscal=0,regions=None,rd=None):
        '''
        Description:
            This fuction is designed to call nucleosome positions
            
        Parameter:
            ofile: a path to the file used to save the positions.
            width: minimal width of positions
            distance: minimal distance between positions, neighboring positions with distance shorter than this value will be merged.
            edge: set to 1 if need to search for position edges, else set to 0
            pcut: a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            fill_gap: fill the gap between two neighboring nuclesomes with a new nucleosome if the gap size is reasonable.
            fill_value: the default value to be set to a filled nucleosome
            calculate_P_value: calculate P value for each position if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
            poscal: set to 1 if need to calculate nucleosome positioning score and P value,else set to 0
        
        Value:
            None.
            
        '''
        outf=open(ofile,mode)
        if mode=='w':
            if poscal>0:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\tfuzziness_score\tfuzziness_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tfuzziness_score\n")
            else:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\n")
        twig=self
        ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log.p=TRUE)/log(10))}''')
        m=self.mean()
        if height==0 and pcut!=0:
            height=int(m+0.5)
            while( (0-float(str(ppois(height,m)).split()[-1])) < pcut):height+=1
        #if not suppress:
        print 'whole genome aveage value is '+str(m)+', use calling cutoff',height
            
        print 'calling summits ...'
        smts=twig.callSummits(width=width,pcut=1,height=height,regions=regions)
        print 'merging summits ...'
        smts.merge(wg=twig,distance=distance)
        if fill_gap:
            print 'filling gaps ...'
            smts.fillgap(wg=twig,height=height,distance=distance)
        
        
        if poscal>0:
            print 'calculating positioning score'
            smts.positioning(twig,rd=rd)
        if edge:print 'searching position edges ...'
        else:print 'saving positions ...'
        width/=2*twig.step
        dic=smts.data
        rhalfdis=distance/2
        halfdis=distance/(2*twig.step)
        for chr in dic:
            lth=twig.chrSize(chr)/twig.step-1
            print chr
            if lth==0:continue
            poses=dic[chr]['p']
            valus=dic[chr]['v']
            if poscal>0:
                positioning=dic[chr]['s']
                ppos=dic[chr]['ppos']
            tlen=poses.size
            if calculate_P_value:
                pvs=ppois(FloatVector(dic[chr]['v']),m)
            i=0
            while i<tlen:
                pos=poses[i]
                start,end=0,0
                if edge==0:start,end=pos-74+(74%self.step)+1,pos+74-(74%self.step)+1
                else:
                    ppp=pos/twig.step
                    p=ppp-1
                    while(start==0):
                        if p<=width:start=1
                        if ppp-p>=halfdis/2:
                            if twig.data[chr][p]==min(twig.data[chr][(p-width):(p+width+1)]):start=p*twig.step+1 
                            elif twig.data[chr][p]>valus[i]:start=p*twig.step+1
                            elif twig.data[chr][p]<height:start=p*twig.step+1
                        elif i>0 and p*twig.step==poses[i-1]:start=p*twig.step+1 
                        p-=1
                    p=ppp+1
                    while(end==0):
                        if (p+width)>=lth:end=lth*twig.step+1
                        elif p-ppp>=halfdis/2:
                            if twig.data[chr][p]==min(twig.data[chr][(p-width):(p+width+1)]):end=p*twig.step+1
                            elif twig.data[chr][p]>valus[i]:end=p*twig.step+1
                            elif twig.data[chr][p]<height:end=p*twig.step+1
                        elif i<(tlen-1) and p*twig.step==poses[i+1]:end=p*twig.step+1 
                        p+=1
                if end>start:
                    if poscal>0:
                        if calculate_P_value:
                            outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\t"+str(positioning[i])+"\t"+str(ppos[i])+"\n")
                        else:outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(positioning[i])+"\n")
                    else:
                        if calculate_P_value:
                            outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\n")
                        else:outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\n")
                i+=1
        outf.close()
        return smts
        

    def fillPositions(self,dic,file,width=40,distance=165,edge=1,pcut=1e-5,height=5,calculate_P_value=1,mode='w',title_line=1,poscal=0,rd=None):
        '''
        Description:
            This fuction is designed to call nucleosome positions
            
        Parameter:
            ofile: a path to the file used to save the positions.
            width: minimal width of positions
            distance: minimal distance between positions, neighboring positions with distance shorter than this value will be merged.
            edge: set to 1 if need to search for position edges, else set to 0
            pcut: a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            fill_gap: fill the gap between two neighboring nuclesomes with a new nucleosome if the gap size is reasonable.
            fill_value: the default value to be set to a filled nucleosome
            calculate_P_value: calculate P value for each position if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
            poscal: set to 1 if need to calculate nucleosome positioning score and P value,else set to 0
        
        Value:
            None.
            
        '''
        twig=self
        smts=Summits()
        smts.data=deepcopy(dic)
        '''
        smts=Summits()
        for line in open(file).readlines()[1:]:
            col=line.split()
            if not smts.data.has_key(col[0]):
                smts.data[col[0]]={}
                smts.data[col[0]]['p']=[]
            smts.data[col[0]]['p'].append(int(col[3])-1)
        for cr in smts.data:
            smts.data[cr]['p'].sort()
            smts.data[cr]['p']=numpy.array(smts.data[cr]['p'])
        '''
        
        smts.fetchValueFromWig(twig)
        if poscal>0:
            print 'calculating positioning score'
            smts.positioning(twig,rd=rd)
            
        outf=open(file,mode)
        
        if mode=='w':
            if poscal>0:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\tfuzziness_score\tfuzziness_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tfuzziness_score\n")
            else:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\n")
        ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log.p=TRUE)/log(10))}''')
        m=self.mean()
        if height==0 and pcut!=0:
            height=int(m+0.5)
            while( (0-float(str(ppois(height,m)).split()[-1])) < pcut):height+=1
        #print 'calling summits ...'
        #smts=twig.callSummits(width=width,pcut=1,height=height,regions=regions)
        #print 'merging summits ...'
        #smts.merge(wg=twig,distance=distance)
        #if fill_gap:
        #    print 'filling gaps ...'
        #    smts.fillgap(wg=twig,height=height,distance=distance)
        if edge:print 'searching position edges ...'
        else:print 'saving positions ...'
        width/=2*twig.step
        dic=smts.data
        rhalfdis=distance/2
        halfdis=distance/(2*twig.step)
        for chr in dic:
            lth=twig.chrSize(chr)/twig.step-1
            print chr
            if lth==0:continue
            poses=dic[chr]['p']
            valus=dic[chr]['v']
            if poscal>0:
                positioning=dic[chr]['s']
                ppos=dic[chr]['ppos']
            tlen=poses.size
            if calculate_P_value:
                pvs=ppois(FloatVector(dic[chr]['v']),m)
            i=0
            while i<tlen:
                pos=poses[i]
                '''
                if valus[i]<height:
                    i+=1
                    continue
                '''
                start,end=0,0
                if edge==0:start,end=pos-74+(74%self.step)+1,pos+74-(74%self.step)+1
                else:
                    ppp=pos/twig.step
                    p=ppp-1
                    while(start==0):
                        if p<=width:start=1
                        if ppp-p>=halfdis/2:
                            if twig.data[chr][p]==min(twig.data[chr][(p-width):(p+width+1)]):start=p*twig.step+1 
                            elif twig.data[chr][p]>valus[i]:start=p*twig.step+1
                            elif twig.data[chr][p]<height:start=p*twig.step+1
                        elif i>0 and p*twig.step==poses[i-1]:start=p*twig.step+1 
                        p-=1
                    p=ppp+1
                    while(end==0):
                        if (p+width)>=lth:end=lth*twig.step+1
                        elif p-ppp>=halfdis/2:
                            if twig.data[chr][p]==min(twig.data[chr][(p-width):(p+width+1)]):end=p*twig.step+1
                            elif twig.data[chr][p]>valus[i]:end=p*twig.step+1
                            elif twig.data[chr][p]<height:end=p*twig.step+1
                        elif i<(tlen-1) and p*twig.step==poses[i+1]:end=p*twig.step+1 
                        p+=1
                if end>start:
                    if poscal>0:
                        if calculate_P_value:
                            outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\t"+str(positioning[i])+"\t"+str(ppos[i])+"\n")
                        else:outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(positioning[i])+"\n")
                    else:
                        if calculate_P_value:
                            outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\n")
                        else:outf.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\n")
                i+=1
        outf.close()
        return smts
        

    def callSummits(self,width=40,pcut=1,height=5,regions=None):
        '''
        Description:
            call occupancy summits using a sliding window.
        
        Parameter:
            width: the width of the sliding window used to call summits.
            pcut:  a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            regions: A set of regions in which the summits are to be defined, regions[chromatin_name][start_position]=end_position
        '''
        if pcut!=1:
            qp=r('''function(p,avg){return(qpois(log(p),avg,lower.tail=FALSE,log=TRUE))}''')
            m=self.mean()
            height=int(str(qp(pcut,m)).split()[-1])
            print 'set summit calling cutoff to',height
        smts=Summits()
        width/=2
        
        if regions==None:
            regions={}
            for cr in self.data:
                regions[cr]={}
                regions[cr][width]=self.data[cr].size*self.step-1
        
        for cr in regions:#self.data:
            print cr
            poses=numpy.array([0])
            valus=numpy.array([0.0])
            dic={}
            lst=self.data[cr]
            step=self.step
            width=width/step
            backwidth=width+1
            #sz=lst.size
            region_poses=regions[cr].keys()
            region_poses.sort()
            num=0
            for p in region_poses:
                i,dlth=p/step,regions[cr][p]/step
                while i<dlth:#lst.size:
                    v=lst[(i-width):(i+backwidth)].max()
                    if v<height:
                        i+=backwidth#continue
                    elif v==lst[i]:
                        #dic[i*step]=lst[i]
                        if num>=poses.size-1:
                            poses.resize(num+1000,refcheck=0)
                            valus.resize(num+1000,refcheck=0)
                        ti,tstart,tend,v1,v2=i,i,i,v,v
                        
                        while v1==v:
                            tstart-=1
                            if tstart>=0:v1=lst[tstart]
                            else:v1,tstart=v-1,0
                        while v2==v:
                            tend+=1
                            if tend<dlth:v2=lst[tend]
                            else:v2,tend=v-1,dlth-1
                        
                        ti=(tstart+tend)/2
                        i=tend
                        poses[num]=ti*step
                        valus[num]=v
                        num+=1
                        i+=backwidth
                        if i<tend:i=tend
                    else:i+=1
            dic['p']=poses[:num]
            dic['v']=valus[:num]
            smts.data[cr]=dic#self[cr].summit(width=width,height=height)#dic
        return smts
    def changeStep(self,step):
        '''
        Description:
            change the step size.
            Note: The new value for each step is determined by sampling an old value within the step!
        
        Parameter:
            step: the new step that is to be setted as.
        
        Value:
            None
        '''
        if step==self.step:return True
        print 'change wiggle step from',self.step,'to',step
        for chr in self.data:
            lth=len(self.data[chr])*self.step
            lst=numpy.array([0.0])
            lst.resize(lth/step+1,refcheck=0)
            for pos in range(0,lth,step):
                try:lst[pos/step]=self.data[chr][pos/self.step]
                except:print lth,pos
            self.data[chr]=lst
        self.step=step
        return True

    def chisqTest(self,wig2,tchr=''):
        '''
        Description:
            do Chi-square test to calculate differential signial for each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance to be compared to
            tchr: specify a chrosome that is to be compared, leave it to '' if want to do for all chrosomes.
            
        Value:
            A Wig class instance that cantain data for the differential signial
        
        '''
        r.options(warn=-2)
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:wig2.changeStep(wig1.step)
        out=Wig(step=wig1.step)
        sum1=wig1.sum()
        sum2=wig2.sum()
        tlen=0
        for chr in wig1.data:tlen+=len(wig1.data[chr])
        donenum=0
        ctime=time()
        for chr in wig1.data:
            if tchr!='' and chr!=tchr:continue
            if wig2.data.has_key(chr):
                wig1.data[chr]+=1
                wig2.data[chr]+=1
                out.data[chr]=numpy.array(0.0)
                lth1=len(wig1.data[chr])
                lth2=len(wig2.data[chr])
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chr].resize(lth,refcheck=0)
                elif lth>lth1:wig1.data[chr].resize(lth,refcheck=0)
                out.data[chr].resize(lth,refcheck=0)
                for p in range(0,lth):
                    #if wig1.data[chr][p]<1:wig1.data[chr][p]=1
                    #if wig2.data[chr][p]<1:wig2.data[chr][p]=1
                    vec=r.c(wig1.data[chr][p],wig2.data[chr][p],sum1-wig1.data[chr][p],sum2-wig2.data[chr][p])
                    if (wig1.data[chr][p]==wig2.data[chr][p]):out.data[chr][p]=0
                    else:
                        test=r['chisq.test'](r.matrix(vec,nrow=2))
                        pvl=float(str(test[2]).split()[-1])
                        if pvl<1e-323:pvl=1e-323
                        if wig1.data[chr][p]>wig2.data[chr][p]:out.data[chr][p]=-log10(pvl)
                        else:out.data[chr][p]=log10(pvl)
                    donenum+=1
                    if time()-ctime>10:
                        ctime=time()
                        print donenum*100.0/tlen,"percent of",tlen,"done"
        wig1.clearEmptyEnd()
        wig2.clearEmptyEnd()
        return out

    def chrSize(self,chr):
        '''
        Description:
            retrive chrosome size by name
            
        Parameter:
            chr: the name of chrosome whose size is to be retrived
        Value:
            Interger value (step size has been multiplied)
            
        '''
        if self.data.has_key(chr):return len(self.data[chr])*self.step
        else:return 0
    def chrSum(self,chr):
        '''
        Description:
            Retrieve the sum of occupancy by chrosome
        Parameter:
            chr: the name of chrosome whose occupancy sum is to be retrieved
        Value:
            Float value
        '''
        return self.data[chr].sum()*self.step
    def clearEmptyEnd(self):
        '''
        Description:
            Clear the 0 values at the end of each chrosome
        Parameter:
            None
        Value:
            None
        '''
        for chr in self.data:
            size=len(self.data[chr])
            while self.data[chr][size-1]==0 and size>0:size-=1#self.data[chr].pop()
            self.data[chr].resize(size,refcheck=0)
    def divideAndLog2(self,wig2):
        '''
        Description:
            divid by wig2 at each data point and then transform the resultant value by log2
        
        Parameter:
            wig2: the Wig class instance that will be used to devide
            
        Value:
            A Wig class instance
        '''
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:
            wig2.changeStep(wig1.step)
        chrs={}
        for chr in wig1.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            lth1=len(wig1.data[chr])
            lth2=len(wig2.data[chr])
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            else:wig1.data[chr].resize(lth2,refcheck=0)
            lth=max(lth1,lth2)
            wig1.data[chr]+=1
            wig2.data[chr]+=1
            wig1.data[chr]/=wig2.data[chr]
            wig1.data[chr][0:lth]=r.log2(FloatVector(wig1.data[chr]))[0:lth]
        return wig1
    def divide(self,wig2):
        '''
        Description:
            divid by wig2 at each data point and then transform the resultant value by log2
        
        Parameter:
            wig2: the Wig class instance that will be used to devide
            
        Value:
            A Wig class instance
        '''
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:
            wig2.changeStep(wig1.step)
        chrs={}
        for chr in wig1.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            lth1=len(wig1.data[chr])
            lth2=len(wig2.data[chr])
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            else:wig1.data[chr].resize(lth2,refcheck=0)
            lth=max(lth1,lth2)
            #wig1.data[chr]+=1
            wig2.data[chr]+=1
            wig1.data[chr]/=wig2.data[chr]
            #wig1.data[chr][0:lth]=r.log2(FloatVector(wig1.data[chr]))[0:lth]
        return wig1
    def dfTest(self,cwig,test='C'):
        '''
        Description:
            Do differential test wit cwig
        Parameter:
            cwig: the Wig class instance which is to be tesed to.
            test: the statistical method that will be used to do the differential test
        Value:
            a Wig class instance containing the differential signal data
        '''
        if test=='C':
            print 'Chi-square test'
            pwig=self.chisqTest(cwig)
            return pwig
        elif test=='F':
            pwig=self.divideAndLog2(cwig)
            return pwig
        elif test=='P':
            print 'Poisson test'
            pwig=self.ppois(cwig)
            return pwig
        elif test=='S':
            print 'subtraction'
            pwig=deepcopy(self)
            pwig.subtract(cwig)
            return pwig
        elif test=='N':
            print "No test method appointed, will not do any differential test.\n"
            return Wig(step=self.step)
        else:
            print "Normalization method "+str(test)+" not applicable now"
            return Wig(step=self.step)

    def fisherTest(self,wig2):
        '''
        Description:
            Do Fisher's exact test to wig2.
            
        Parameter:
            wig2: the Wig class instance which is to be tesed to.
            
        Value:
            a Wig class instance containing the differential signal data
        '''
        r.options(warn=-1)
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        out=Wig(step=wig1.step)
        sum1=int(self.sum())
        sum2=int(wig2.sum())
        tlen=0
        for chr in self.data:tlen+=len(self.data[chr])
        donenum=0
        for chr in self.data:
            if wig2.data.has_key(chr):
                out.data[chr]=numpy.array(0.0)
                lth1=len(self.data[chr])
                lth2=len(wig2.data[chr])
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chr].resize(lth1,refcheck=0)
                elif lth>lth1:self.data[chr].resize(lth2,refcheck=0)
                out.data[chr].resize(lth,refcheck=0)
                for p in range(0,lth):
                    vec=r.c(int(self.data[chr][p]),int(wig2.data[chr][p]),sum1-int(self.data[chr][p]),sum2-int(wig2.data[chr][p]))
                    if (self.data[chr][p]==wig2.data[chr][p]):out.data[chr][p]=1
                    else:
                        test=r['fisher.test'](r.matrix(vec,nrow=2))
                        pvl=float(str(test).split()[17])
                        if pvl<1e-323:pvl=1e-323
                        if self.data[chr][p]>wig2.data[chr][p]:out.data[chr][p]=-log10(pvl)
                        else:out.data[chr][p]=log10(pvl)
                    donenum+=1
                    if donenum%1000==0:print donenum*100.0/tlen,"percent of",tlen,"done"
        self.clearEmptyEnd()
        wig2.clearEmptyEnd()
        return out
    def foldChange(self,fold):
        '''
        Description:
            Do fold change at each data point.
            
        Parameter:
            fold: the value that will be multiplied by each data point
            
        Value:
            None
        '''
        for chr in self.data:self.data[chr]*=fold

    def getChrs(self):
        '''
        Description:
            Retrive all chrosome names
        Parameter:
            None
        Value:
            a list of chrosome names
        '''
        return self.data.keys()
        
    def gsize(self):
        '''
        Description:
            Calculate genome size
        Parameter:
            None
        Value:
            Interger value
        '''
        lth=0
        for chr in self.data:
            lth+=self.chrSize(chr)#self.step*self.data[chr].size
        return lth
    def histogram(self,bnum=100000,nonzero_end=False):
        ma,mi=self.maxmin(nonzero=nonzero_end)
        counts,bins=numpy.histogram([],bins=bnum,range=(mi,ma))
        for cr in self.data:counts+=numpy.histogram(self.data[cr],bins=bnum,range=(mi,ma))[0]
        return [counts,bins]

    def load(self,file,gfile,suppress=False):
        '''
        Description:
            Load Wig class instance from Wiggle format file
        Parameter:
            file: a path to the file containing the data
        Value:
            None
        '''
        if file[-3:]=='wig':fi=open(file)
        if file[-6:]=='wig.gz':
            import gzip
            fi=gzip.open(file)
        for line in fi:
            if line[0]=='t':continue
            elif line[0:9]=='fixedStep':
                self.loadFixed(file,gfile=gfile,suppress=suppress)
                return
            elif line[0:3]=='var':
                self.loadVar(file,gfile=gfile,suppress=suppress)
                return
            else:print "Load failure: format not recoganized!",line
    
    def loadFixed(self,file,gfile,suppress=False):#add by Kaifu April 4, 2012
        '''
        Description:
            load data from Fixed wiggle format file
        Parameter:
            file:a path to the file containging the data
            suppress: suppress waring message? True or False
        Value:
            None
        '''
        sss=time()
        
        ########## ---start---add by kaifu on Aug 14, 2012 ##########
        if self.step<1:
            tempf=open(file)
            if not suppress: print 'detecting step size ...'
            while self.step<1:
                line=tempf.readline()
                if line[0]=='f':
                    col=line.split()
                    for term in col[1:]:
                        kv=term.split("=")
                        if kv[0]=='step':self.step=int(kv[1])
        ########## ---end---add by kaifu on Aug 14, 2012 ##########
        
        if not suppress: print 'parsing from',file
        chr=''
        chrs=['chr1']############ test #############
        pn=1
        pos=-1
        for line in open(file):
            if line[0]=='t':continue
            elif line[0]=='f':
                col=line.split()
                for term in col[1:]:
                    kv=term.split("=")
                    if kv[0]=='chrom':
                        newchr=kv[1]
                    elif kv[0]=='start':newstart=int(kv[1])
                    elif kv[0]=='step':instep=int(kv[1])
                
                if not suppress:print newchr,'start from position ',newstart#/instep
                if chr != '':
                    if not self.data.has_key(chr):self.data[chr]=numpy.array([0.0])
                    nplst=numpy.array(lst)#[0:]
                    if instep!=self.step:
                        lth=len(lst)*(instep*1.0/self.step)
                        nplst=numpy.array([0.0])
                        nplst.resize(int(lth),refcheck=0)
                        for pos in range(len(lst)):
                            #nplst[pos*instep/self.step]+=lst[pos]########## deleted by kaifu on Sep 6, 2012 ##########
                            if instep<self.step:nplst[pos*instep/self.step]+=lst[pos] ########## add by kaifu on Sep 6, 2012 ##########
                            else:
                                for tpos in range(pos*instep/self.step,(pos+1)*instep/self.step):nplst[tpos/self.step]+=lst[pos]*1.0*self.step/instep ########## add by kaifu on Sep 6, 2012 ##########
                    self.data[chr].resize(start/self.step+len(nplst),refcheck=0)
                    self.data[chr][(start/self.step):]=nplst
                lst=[]
                chr=newchr
                start=newstart-1
            else:
                try:
                    #if line=='nan\n':line='0\n'
                    lst.append(float(line.split()[0]))
                except:
                    if line[:3]=='nan':lst.append(0.0)
                    else:print 'wrong line:',line[:-1]
                    continue
                #print line.split()[0],lst[-1]
                #lst.append(value)
        #if chr != '':
        if not self.data.has_key(chr):self.data[chr]=numpy.array([0.0])
        nplst=numpy.array(lst)#[0:]
        if instep!=self.step:
            lth=len(lst)*(instep*1.0/self.step)
            nplst=numpy.array([0.0])
            nplst.resize(int(lth),refcheck=0)
            for pos in range(len(lst)):
                #nplst[pos*instep/self.step]=lst[pos] ########## deleted by kaifu on Sep 6, 2012 ##########
                if instep<self.step:nplst[pos*instep/self.step]+=lst[pos] ########## add by kaifu on Sep 6, 2012 ##########
                else:
                    for tpos in range(pos*instep/self.step,(pos+1)*instep/self.step):nplst[tpos/self.step]+=lst[pos]*1.0*self.step/instep ########## add by kaifu on Sep 6, 2012 ##########
        self.data[chr].resize(start/self.step+len(nplst),refcheck=0)
        self.data[chr][(start/self.step):]=nplst

    def loadVar(self,file,gfile,suppress=False):
        '''
        Description:
            load data from Fixed wiggle format file
        Parameter:
            file:a path to the file containging the data
        Value:
            None
        '''
        
        ########## ---start---add by kaifu on Aug 14, 2012 ##########
        if self.step<1:
            print 'set step size to 10'
            self.step=10
        ########## ---end---add by kaifu on Aug 14, 2012 ##########
        print 'parsing from',file
        step=self.step

        starttime=time()
        chr,size='',0
        if file[-3:]=='wig':fi=open(file)
        if file[-6:]=='wig.gz':
            import gzip
            fi=gzip.open(file)
        for line in fi:
            #print line
            if line[0]=='t':continue
            elif line[0]=='v':
                #if chr != '':
                #    if not self.data.has_key(chr):self.data[chr]=numpy.array([0.0])#lst
                #    print line[:-1]
                if not suppress: print line[:-1]
                col=line.split()
                right=0
                for term in col[1:]:
                    kv=term.split("=")
                    if kv[0]=='chrom':
                        right+=1
                        chr=kv[1]
                        #print chr,#line[:-1]
                        if not self.data.has_key(chr):
                            self.data[chr]=numpy.array([0.0])#lst
                        lst=self.data[chr]#chr=kv[1]
                        size=lst.size
                    elif kv[0]=='span':
                        inspan=int(kv[1])
                        right+=1
                if right<2:
                    print 'wrong format:',line,'chrom and span must be provided'
            else:
                col=line.split()
                tstart,value=int(col[0]),float(col[1])
                tend=(tstart+inspan)/step
                tstart=tstart/step
                if tend>=size:
                    #while(tend>=size):
                    size=tend+1000
                    lst.resize(size,refcheck=0)
                for pos in range(tstart,tend):lst[pos]=value
        self.clearEmptyEnd()
        #print 'time cost:',time()-starttime
    def maxmin(self,nonzero=False):
        ma,mi=0,0
        for cr in self.data:
            tma,tmi=self.data[cr].max(),self.data[cr].min()
            if nonzero:
                tpos=self.data[cr].nonzero()
                tma,tmi=self.data[cr][tpos].max(),self.data[cr][tpos].min()
                if mi==0:mi=tmi
                if ma==0:ma=tma
            if tma>ma:ma=tma
            if tmi<mi:mi=tmi
        return [ma,mi]
    
    def mean(self):
        '''
        Description:
            Calculate the mean occupancy value
        Parameter:
            None
        Value:
            float value
        '''
        return self.sum()*1.0/self.gsize()
    def multiply(self,wig2):
        '''
        Description:
            multiply by wig2 at each data point
        Parameter:
            wig2: the Wiggle class instance to multiply by
        Value:
            an Wiggle class instance
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        wg=Wig(step=self.step)
        chrs={}
        for chr in self.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            if not self.data.has_key(chr):self.addChr(chr)
            if not wig2.data.has_key(chr):wig2.addChr(chr)
            lth1=self.data[chr].size
            lth2=wig2.data[chr].size
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            elif lth1<lth2:self.data[chr].resize(lth2,refcheck=0)
            wg.data[chr]=self.data[chr]*wig2.data[chr]
        return wg
    def non0to1(self):
        '''
        change all non-zero value to value 1
        #can be further improved by using the numpy.nonzero() function
        '''
        for cr in self.data:
            a=self.data[cr]
            a[a.nonzero()]=1

    def non0size(self):
        '''
        the number of data point with a non-zero value.
        '''
        size=0
        for cr in self.data:
            size+=self.data[cr].nonzero()[0].size
        return size*self.step

    def percentile(self, p=[0,25,50,75,100],bnum=100000,nonzero_end=False):
        '''
        bnum: number of histogram bins to be calculated between maximal value and minimal value
        '''
        hs=self.histogram(bnum=bnum,nonzero_end=nonzero_end)
        total=hs[0].sum()
        p=deepcopy(p)
        p.sort()
        plth=len(p)
        i,j=0,0
        count=0
        while j<bnum:
            count+=hs[0][j]
            while count*100.0/total>=p[i] and i<plth:
                p[i]=hs[1][j]*(100-p[i])/100.0+hs[1][j+1]*(p[i])/100.0
                i+=1
                if i>=plth:
                    p.append(110)
                    j=bnum
            j+=1
        p=p[:plth]
        #print p
        return p
    
    def pop(self,k):
        '''
        Description:
            remove a chrosome 
        Parameter:
            k: the name of the chrosome that is to be removed
        Value:
            None
        '''
        return self.data.pop(k)
    def power(self,p):
        '''
        Description:
            Self multiply by p times
        Parameter:
            p: the times to do self multiply
        Value:
            None
        '''
        chrs={}
        for chr in self.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            self.data[chr]=self.data[chr]**p
        return True

    def ppois(self,wig2,tchr=''):
        '''
        Description:
            do Poisson test to calculate differential signial for each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance to be compared to
            tchr: specify a chrosome that is to be compared, leave it to '' if want to do for all chrosomes.
            
        Value:
            A Wig class instance that cantain data for the differential signial
        
        '''
        pp=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE)/log(10))}''')###### add by kaifu on Oct 18,2011
        wig1=self#deepcopy(self)
        if wig1.step!=wig2.step:wig2.changeStep(wig1.step)
        out=Wig(step=wig1.step)
        donenum=0
        for chr in wig1.data:
            print chr
            if tchr!='' and chr!=tchr:continue
            if wig2.data.has_key(chr):
                out.data[chr]=numpy.array(0.0)
                lth1=wig1.data[chr].size
                lth2=wig2.data[chr].size
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chr].resize(lth,refcheck=0)
                elif lth>lth1:wig1.data[chr].resize(lth,refcheck=0)
                out.data[chr].resize(lth,refcheck=0)
                tstr1=numpy.array(0.0)
                tstr2=numpy.array(0.0)
                tstr1.resize(lth,refcheck=0)
                tstr2.resize(lth,refcheck=0)
                p=0
                while p<lth:
                    if wig1.data[chr][p]>wig2.data[chr][p]:
                        if wig1.data[chr][p]<1:tstr1[p]=1
                        else:tstr1[p]=wig1.data[chr][p]
                        if wig2.data[chr][p]<1:tstr2[p]=1
                        else:tstr2[p]=wig2.data[chr][p]
                    else:
                        if wig1.data[chr][p]<1:tstr2[p]=1
                        else:tstr2[p]=wig1.data[chr][p]
                        if wig2.data[chr][p]<1:tstr1[p]=1
                        else:tstr1[p]=wig2.data[chr][p]
                    p+=1
                frglen=1000000   # Fragmentation was implemented by Kaifu Chen on Feb 22, 2013. The reason is that R can not allocate vector space larger than 1Gb.
                frags={}
                for end in range(frglen,lth,frglen):frags[end-frglen]=end
                frags[lth-lth%frglen]=lth
                starts=frags.keys()
                starts.sort()
                for start in starts :
                    end=frags[start]
                    print 'region:',start,end
                    result=pp(FloatVector(tstr1[start:end]),FloatVector(tstr2[start:end]))
                    p=start
                    while p<end:
                        if tstr1[p]==tstr2[p]:out.data[chr][p]=0
                        elif wig1.data[chr][p]>wig2.data[chr][p]:out.data[chr][p]=0-result[p-start]
                        elif wig1.data[chr][p]<wig2.data[chr][p]:out.data[chr][p]=result[p-start]
                        p+=1
        return out

    def regionWithinValueRange(self, lowValue=None,highValue=None):
        '''
        Each value between lowValue and highValue will be set to 1, other value will be set to 0
        #may be further improved by using the numpy.where() or numpy.nonzero() function
        '''
        twig=deepcopy(self)
        if lowValue==None and highValue==None:
            twig.foldChange(0)
            return twig
        ma,mi=twig.maxmin()
        if lowValue==None:lowValue=mi-1
        if highValue==None:highValue=ma+1
        for cr in twig.data:
            twig.data[cr][numpy.where([((self.data[cr]<lowValue) | (self.data[cr]>highValue))])[1]]=0        
            twig.data[cr][numpy.where([((self.data[cr]>=lowValue) & (self.data[cr]<=highValue))])[1]]=1       
        return twig
    
    
    
    def resizeChr(self,chr,size):
        '''
        Description:
            change chrosome size
        Parameter:
            chr: the name of the chrosome whose size is going to be changed
            size: the new size to be setted for the chrosome
        Value:
            None
        '''
        self.data[chr].resize(size/self.step,refcheck=0)
    def rvNeg(self):
        '''
        Description:
            set negative value at each data point to 0
        Parameter:
            None
        Value:
            None
        '''
        for chr in self.data:
            end=self.chrSize(chr)/self.step
            t=self.data[chr]
            self.data[chr]=((t**2)**0.5+t)/2


    def save(self,file,format="fixed",step=None,suppress=False):
        '''
        Description:
            Save data to wiggle format file
        Parameter:
            file: a path to the output file
            format: the format of the output wiggle file, could be 'fixed' or 'var'
            step: the step size of the ouput wiggle file
        Value:
            None
        '''
        if step==None:step=self.step
        if not suppress:print 'saving  to '+file
        outf=open(file,"w")
        if format==None:#add by Kaifu on Feb 10, 2014
            from random import randint
            chrs=self.data.keys()
            lth=len(chrs)-1
            count,total=0,1000
            for i in range(total):
                chr=chrs[randint(0,lth)]
                val=self.data[chr][randint(0,sef.data[chr].size-1)]
                if val!=0:count+=1
            if count*1.0/total>=0.5:format='fixed'
            else:format='var'
            print 'Will save in .'+format+' format'
        if format=='wiq':
            tlen=0
            for cr in self.data:tlen+=self.data[cr].size
            from numpy.random import seed,randint
        acculen=0
        for chr in self.data:
            #if self.data[chr].sum()==0:continue
            if not suppress:print chr
            if format=='fixed':
                outf.write("fixedStep chrom="+chr+" start=1  step="+str(step)+" span="+str(step)+"\n")
                lth=len(self.data[chr])
                ot=[]
                for v in self.data[chr]:ot.append(str(v))#str(self.data[chr][i]))
                outf.write('\n'.join(ot)+"\n")
            elif format=='var':
                lth=len(self.data[chr])
                outf.write("variableStep chrom="+chr+" span="+str(step)+"\n")
                #for i in range(0,lth):
                i=0
                while i<lth:
                    if self.data[chr][i]!=0:outf.write(str(i*step)+'\t'+str(self.data[chr][i])+'\n')
                    i+=1
            elif format=='wiq':
                lth=len(self.data[chr])
                #outf.write("variableStep chrom="+chr+" span="+str(step)+"\n")
                #for i in range(0,lth):
                i=0
                while i<lth:
                    #if self.data[chr][i]!=0:outf.write(str(i*step)+'\t'+str(self.data[chr][i])+'\n')
                    acculen+=1
                    seed(acculen)
                    outf.write(str(self.data[chr][i])+'\t'+str(randint(0,tlen))+'\t'+chr+'\t'+str(i*step)+'\n')
                    i+=1
            else:
                print ',format not recogonized, will be saved in fixed format,',
                outf.write("fixedStep chrom="+chr+" start=1  step="+str(step)+" span="+str(step)+"\n")
                lth=len(self.data[chr])
                for i in range(0,lth):outf.write(str(self.data[chr][i])+"\n")
        outf.close()
        print 'completed'
    def saveChr(self,file,chr=None,lth=None,format="fixed",step=10):
        '''
        Description:
            Save data to wiggle format file by chrosome name
        Parameter:
            file: a path to the output file
            format: the format of the output wiggle file, could be 'fixed' or 'var'
            step: the step size of the ouput wiggle file
            chr: the name of the chrosome that is going to be saved
            lth: the length of the chrosome to be saved, start from 1 to lth
        Value:
            None
        '''
        tchr=chr
        print 'saving wig to',file,
        outf=open(file,"w")
        for chr in self.data:
            if chr!=tchr:continue
            if format=='fixed':
                outf.write("fixedStep chrom="+chr+" start=1  step="+str(step)+" span="+str(step)+"\n")
                if lth==None:lth=len(self.data[chr])*self.step
                for i in range(0,lth,step):outf.write(str(self.data[chr][i/self.step])+"\n")
        outf.close()
        print 'completed'

    def sd(self):
        '''
        Description:
            Calculate standard deviation of occupancy
        Parameter:
            None
        Value:
            None
        '''
        sz=self.size()
        avg=self.sum()/sz
        sqm=0
        for chr in self.data:
            for v in self.data[chr]:sqm+=(v-avg)*(v-avg)
        sqm/=sz
        return sqrt(sqm)
    

    def sizeAdjust(self,gfile):
        '''
        Description:
            Adjust the size of each chrosome.
        
        Parameter:
            gfile: path to the file containing the size of each chrosome, each line in the file would be in the format "chrosome_name size", in which size is an integer value, and chrosome_name should contain no empty space
        
        value:
            None.
        
        '''
        sizes={}
        for line in open(gfile):
            col=line.split()
            sizes[col[0]]=int(col[1])
        for chr in self.data:
            if not sizes.has_key(chr):self.data.pop(chr)
            else:
                print self.data[chr].size,
                self.data[chr].resize(sizes[chr]/self.step,refcheck=0)
                print self.data[chr].size
    def smooth(self,lmd=100):
        '''
        Description:
            Smooth occupancy by a sliding window
        Parameter:
            the size of the smooth window
        Value:
            None
        '''
        print 'smooth width:',lmd
        ss=time()
        lmd=int(lmd/self.step)
        if lmd<=0:return True
        hlmd=lmd/2
        tlmd=lmd-hlmd
        wg2=deepcopy(self.data)
        self.foldChange(0.0)
        wg1=self.data
        for chr in wg1:
            lth=wg1[chr].size
            if tlmd!=0-hlmd:
                for p in range(0-hlmd,tlmd):wg1[chr][hlmd:(lth-tlmd)]+=wg2[chr][(hlmd+p):(lth-tlmd+p)]
            else:
                for p in range(0-hlmd+1,tlmd):wg1[chr][hlmd-1:(lth-tlmd)]+=wg2[chr][(hlmd-1+p):(lth-tlmd+p)]
                for p in [0-hlmd,tlmd]:wg1[chr][hlmd:(lth-tlmd)]+=wg2[chr][(hlmd+p):(lth-tlmd+p)]*0.5
            wg1[chr][hlmd:(lth-tlmd)]/=(lmd)*1.0
        self.data=wg1
        return True
    def sqrt(self):
        '''
        Description:
            translate each data point to its square root value
        Parameter:
            None
        Value:
            None
        '''
        chrs={}
        for chr in self.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            self.data[chr]=numpy.sqrt(self.data[chr])
        return True
    def subtract(self,wig2):
        '''
        Description:
            Subtract wig2 at each data point
        Parameter:
            wig2: the Wig class instance to be subtracted
        Value:
            None
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        chrs={}
        for chr in self.data:chrs[chr]=1
        for chr in wig2.data:chrs[chr]=1
        chrs=chrs.keys()
        for chr in chrs:
            lth1=len(self.data[chr])
            lth2=len(wig2.data[chr])
            if lth1>lth2:wig2.data[chr].resize(lth1,refcheck=0)
            else:self.data[chr].resize(lth2,refcheck=0)
            self.data[chr]-=wig2.data[chr]
        wig2.clearEmptyEnd()
        return True
    def sum(self):
        '''
        Description:
            return the sum of occupancy across the whole genome
        Parameter:
            None
        Value:
            None
        '''
        v=0
        for chr in self.data:v+=numpy.sum(self.data[chr])*self.step
        return v

    def sumWithinValueRange(self, lowValue,highValue):
        '''
        add values between lowValue and highValue, ignore other values
        '''
        twig=self#deepcopy(self)
        total=0
        for cr in twig.data:
            size=twig.data[cr].size()
            i=0
            while i<size:
                if twig.data[cr][i]>=lowValue and twig.data[cr][i]<=highValue:total+=twig.data[cr][i]
                i+=1        
        return total
    
    def sam_coverage(self,sam_file,step=None):
        if step==None:step=self.step
        else:self.step=step
        infile=open(sam_file)#os.popen('samtools view -XS '+sam_file)
        line=infile.readline()
        hlines=[]
        while line[0]=='@':
            hlines.append(line)
            col=line.split()
            if col[0]=='@SQ':
                chr,clen=col[1][3:],int(col[2][3:])
                self.data[chr]=numpy.array([0.0])
                self.data[chr].resize(clen/step,refcheck=0)
            line=infile.readline()
            
        infile=open(sam_file)
        for i in range(len(hlines)):infile.readline()
        for line in infile:
            col=line[:-1].split('\t')
            chr=col[2]
            t1=re.findall('\d+',col[5])#.split('\d+'))
            t2=re.findall('\D+',col[5])
            start=int(col[3])
            for i in range(len(t2)):
                if t2[i]=='M' :
                    end=start+int(t1[i])
                    self.data[chr][start/step:end/step]+=1
                elif t2[i]=='D':start=start+int(t1[i])
                elif t2[i]=='N':start=start+int(t1[i])
                
                '''
                elif t2[i]=='S':continue
                elif t2[i]=='H':continue
                elif t2[i]=='P':continue
                elif t2[i]=='I':continue
                '''
    def std(self):
        v=0
        m=self.mean()
        gs=0
        for chr in self.data:
            ts=self.data[chr].size
            i=0
            while i<ts:
                s=self.data[chr][i]-m
                v+=s*s
                i+=1
            gs+=ts
        v=sqrt(v/gs)
        return v
if __name__ == "__main__":
    import sys,re,os
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.

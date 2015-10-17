#!/usr/bin/env python

from copy import deepcopy
from wig import Wig
from wigs import Wigs
import argparse,sys,os
from lib import peakSelectorByValue,peakSelectorByGeneStructure,batchOccAroundPoints,occAroundPoints,plot,batchOccPSD,retrieve_peaks_by_value,batchPeakDistanceDistribution,batchPeakValDistribution
from rpy2.robjects import r


def printHelp():
    print '\ndantools version 0.2.1'
    print '\nFor help information, try:\npython path/dantools.py <command>'
    print 'Command:'
    print '\twigAnalysis:\tanalyze wiggle format data.'
    print '\tpeaksAnalysis:\tstatistics of nucleosome peaks.'
    print '\tpeakSelector:\tselect a subset of peaks.'
    print '\tretrievePeakValuesAtRanks:\treturn peak values at the sepcified ranks.'
    print '\nKaifu Chen, et al. chenkaifu@gmail.com, Li lab, Biostatistics department, Dan L. Duncan cancer center, Baylor College of Medicine.'
    print ''
    
    
def wiggleAnalysis(args,wgs,pos_neg='', outname='',outmode='w'):
    rcode=''
    #rcode='par(mfrow=c('+str(nrow)+','+str(ncol)+'))\n'#+rcode
    gfiles=args.genefile_paths.split(',')
    wigaliases=args.wigfile_aliases.split(',')
    if args.genefile_aliases==None:gfnames=args.genefile_paths.split(',')
    else:gfnames=args.genefile_aliases.split(',')
    sites=args.genomic_sites.split(',')
    if 'TSS' in sites:
        print '\nprofiling for Transcript Start Sites (TSS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.tss',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                else :dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.tss',outmode='a',chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.tss',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.tss.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.tss',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    if 'TTS' in sites:
        print '\nprofiling for Transcription Terminal Sites (TTS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.tts',outmode='a',chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.tts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.tts.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.tts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
        
    if 'CSS' in sites:
        print '\nprofiling for Coding Start Sites (CSS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.css',outmode='a',chrColID=1,nameColID=0,posColIDpos=5,posColIDneg=6,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.css',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.css.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.css',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    
    if 'CTS' in sites:
        print '\nprofiling for Coding Terminal Sites (CTS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.cts',outmode='a',chrColID=1,nameColID=0,posColIDpos=6,posColIDneg=5,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.cts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.cts.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.cts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    
    if 'ESS' in sites:
        print '\nprofiling for Exon Start Sites (ESS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.ess',outmode='a',chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.ess',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.ess.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.ess',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    if 'ETS' in sites:
        print '\nprofiling for Exon Terminal Sites (ETS)'
        d={}
        if len(wgs.keys())>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg+'.ets',outmode='a',chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',step=args.bin_size,lines=glines,heatMap=True,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.ets',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.ets.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.ets',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    return rcode

def runWiggleAnalysis():
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the main function wiggleAnalysis(), or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\npython dantools.py wigAnalysis <wiggle_file_paths> <gene_file_paths> [optional arguments]\n\nfor more help, please try: python dantools.py peaksAnalysis -h\n"
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython dantools.py wigAnalysis <wiggle_file_paths> <gene_file_paths>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    parser.add_argument('command',default=None,\
                        help="This tells dantools the job that is to be done, could be wigAnalysis, peaksAnalysis, peakSelector, or retrievePeakValuesAtRanks")
    parser.add_argument('wigfile_paths',default=None,\
                        help="Paths to wiggle format data sets. \
                        each path directs to a wiggle format file, \
                        paths must be separated by the comma ','. e.g. file_a.wig,file_b.wig,file_c.wig")
    parser.add_argument('genefile_paths',default=None,\
                        help="Paths to files each contains a gene sets. \
                        paths must be separated by the comma ','.\
                        Each gene set file must contain at least the following columns ordered as: name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds, \
                        we suggest to download gene set from the UCSC tables at http://genome.ucsc.edu/cgi-bin/hgTables?command=start ")
    parser.add_argument('--genomic_sites',dest='genomic_sites',metavar='',default='TSS,TTS,CSS,CTS,ESS,ETS',\
                        help="The genomic site names to be analyzed, choose one or several from transcription start site (TSS), \
                        transcription terminal site (TTS), coding start site (CSS), coding terminal site (CTS), exon start site (ESS), \
                        or exon terminal site (ETS). Names must be separated by the comma ',' ")
    parser.add_argument('--periodicity',dest='periodicity',metavar='',type=int,default=1,\
                        help="Set to 1 to do power spectrum density (PSD) analysis, else set to 0 to cancel this analysis. This function checks the strength of periodicity in each wiggle data set.")
    parser.add_argument('--wigfile_aliases',dest='wigfile_aliases',metavar='',default=None,\
                        help="A set of short aliases for the wiggle files. \
                        Aliases must be separated by the comma ',' and arranged in the same order as the wiggle files. The data set in each wiggle file \
                        will be represent by the alias in the output files, or by the file name when no alias is provided")
    parser.add_argument('--genefile_aliases',dest='genefile_aliases',metavar='',default=None,\
                        help="A set of short aliases for the gene files. \
                        Aliases must be separated by the comma ',' and be arranged in the same order as the gene files. Each gene \
                        file will be represented by its alias in the output files, or by the file name when no alias is provided.")
    parser.add_argument('--name', dest="name",metavar='',default='wigAnalysis',\
                        help="A name for the experiment.")
    parser.add_argument('--vcal', dest="vcal",metavar='',default='mean',\
                        help="The method to calculate plot value at each position relative to a genomic site (e.g. TSS) in a gene group, could be 'median' or 'mean'")
    parser.add_argument('--excludeP',dest="excludeP",metavar='',default=1,type=float,\
                        help="Exclude the extremely outgroup (low or high) values by this percentage.")
    parser.add_argument('--pos_neg',dest="pos_neg",metavar='',default=0,type=int,\
                        help="Set to 0, 1, -1, 2, or 3 if want to do analysis for \
                        positive and negative values together, positive values only, \
                        negative values only, positive and negative values seperately, \
                        or all kinds of analysis (together, positive, negative),\
                        this will be useful when we are analyze the differential signal between two samples.")
    parser.add_argument('--bin_size',dest="bin_size",metavar='',default=10,type=int,\
                        help="Bin size to be used for plotting, it is suggested to be the same as the step or span size in the input wiggle files.")
    parser.add_argument('--flank_up',dest="flank_up",metavar='',default=1500,type=int,\
                        help="How far to calculate from the up stream of each genomic site (e.g., TSS).")
    parser.add_argument('--flank_dn',dest="flank_dn",metavar='',default=3000,type=int,\
                        help="How far to calculate to the down stream of each genomic site (e.g., TSS).")
    parser.add_argument('--plot_row',dest="plot_row",metavar='',default=2,type=int,\
                        help="Number of rows to plot on each page.")
    parser.add_argument('--plot_column',dest="plot_column",metavar='',default=2,type=int,\
                        help="Number of columns to polt on each page.")
    parser.add_argument('--plot_xmin',dest="plot_xmin",metavar='',default=None,type=float,\
                        help="Minimal value on the x axis of each plot.")
    parser.add_argument('--plot_xmax',dest="plot_xmax",metavar='',default=None,type=float,\
                        help="Maximal value on the x axis of each plot.")
    parser.add_argument('--plot_ymin',dest="plot_ymin",metavar='',default=None,type=float,\
                        help="Minimal value on the y axis of each plot.")
    parser.add_argument('--plot_ymax',dest="plot_ymax",metavar='',default=None,type=float,\
                        help="Maximal value on the y axis of each plot.")
    parser.add_argument('--plot_xlab', dest="plot_xlab",metavar='',default='Relative distance',\
                        help="The label on the x axis.")
    parser.add_argument('--plot_ylab', dest="plot_ylab",metavar='',default='Average signal value',\
                        help="The label on the y axis.")
    parser.add_argument('--plot_colors', dest="plot_colors",metavar='',default='black,gray,red,blue,orange,purple,skyblue,cyan,green,blue4,darkgoldenrod',\
                        help="The colors to be used in the plot.")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0
    elif len(sys.argv)>=4: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print "\nfor more help, please try: python dantools.py peaksAnalysis -h\n"
            return 0
    else:
        print "\nfor help, please try: python dantools.py peaksAnalysis -h\n"
        return 0 
    
    print '\ncommand:\npython'," ".join(sys.argv) # print the command line, this let the user to keep a log and remember what parameters they specified
    #print args
    print '\n\nparsing wiggle format data ...'
    if args.wigfile_aliases==None:args.wigfile_aliases=args.wigfile_paths
    wgs=Wigs()
    wigpaths=args.wigfile_paths.split(',')
    wigaliases=args.wigfile_aliases.split(',')
    if len(wigpaths)!=len(wigaliases):
        print 'Error: wigfile alias count not equal to path count!\n'
        #parser.print_help()
        print ''
        return
    else:
        for i in range(len(wigpaths)):
            wgs.set(wigaliases[i],Wig(wigpaths[i],suppress=True))
            print wigaliases[i],wgs.data[wigaliases[i]].sum()

    print '\n\nprofiling ...'
    rcode=''
    if args.pos_neg==0 or args.pos_neg==3:
        print '\n\nprofiling  for positive and negative values together...'
        rcode+=wiggleAnalysis(args,wgs,pos_neg='',outname=args.name,outmode='w')
    
    if args.pos_neg==1 or args.pos_neg==2 or args.pos_neg==3:
        print '\n\nprofiling  for positive values ...'
        pwgs=deepcopy(wgs)
        for k in pwgs.keys():pwgs.get(k).rvNeg()
        if args.pos_neg==1 or args.pos_neg==2: rcode+=wiggleAnalysis(args,pwgs,pos_neg='.pos',outname=args.name,outmode='w')
        else:rcode+=wiggleAnalysis(args,pwgs,pos_neg='.pos',outname=args.name,outmode='a')
        
    if args.pos_neg==-1 or args.pos_neg==2 or args.pos_neg==3:
        print '\n\nprofiling  for negative values ...'
        nwgs=deepcopy(wgs)
        for k in nwgs.keys():
            nwgs.get(k).foldChange(-1.0)
            nwgs.get(k).rvNeg()
            nwgs.get(k).foldChange(-1.0)
        if args.pos_neg==2 or args.pos_neg==3: rcode+=wiggleAnalysis(args,nwgs,pos_neg='.neg',outname=args.name,outmode='a')
        else:rcode+=wiggleAnalysis(args,nwgs,pos_neg='.neg',outname=args.name,outmode='w')
    
    
    if args.periodicity!=0:
        print '\n\ncalculating periodicity strength...'
        dic=batchOccPSD(wgs,outname=args.name+'.singalPeriodicity')
        rcode+=plot(dic=dic,outname='',main='singalPeriodicity',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab='Periodicity Length (bp)',ylab='Strength',colors=args.plot_colors.split(','))

    #rcode='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'+rcode
    rcode='pdf("'+args.name+'.pdf")\n'+rcode
    rcode+='dev.off()\n'
    fo=open(args.name+'.R','w')
    fo.write(rcode)
    fo.close()
    r(rcode)
    print ''
    
def runPeaksAnalysis():
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the main functions related to peaks analysis, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\npython dantools.py peaksAnalysis <peakfile_paths> <peakfile_aliases> [optional arguments]\n\nfor more help, please try: python dantools peaksAnalysis -h\n"
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython dantools.py peaksAnalysis <peakfile_paths> <peakfile_aliases>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    parser.add_argument('command',default=None,\
                        help="This tells dantools to run peaksAnalysis")
    parser.add_argument('peakfile_paths',default=None,\
                        help="Paths to the peak files. Paths must be separated by the comma ',', \
                        each peak file should be in the default output format of DANPOS, e.g. the file *.allPeaks.xls or the files pooled/*.peaks.xls. \
                        Each *.allPeaks.xls file contains peaks from a control sample and a treat sample, \
                        each pooled/*.peaks.xls file contains peaks from a single sample.\
                        See DANPOS documentation for more information.\
                        e.g. a-b.allpeaks.xls,c.peaks.xls,d-e.peaks.xls")
    parser.add_argument('peakfile_aliases',default=None,\
                        help="A set of aliases for the peak files. \
                        Each *.allPeaks.xls file would have a pair of aliases in the \
                        format a:b, 'a' represent the sample name of treat peaks, \
                        'b' represents the sample name of control peaks. \
                        Each *.peaks.xls file generated by DANPOS in the directory 'pooled/' would have a single alias.\
                        Alias for different peak files must be separated by the comma ',', \
                        and must be arranged in an order corresponding to the peakfile_paths.\
                        e.g. a:b,c,d:e")
    parser.add_argument('--name', dest="name",metavar='',default='peaksAnalysis',\
                        help="A name for the experiment.")
    
    parser.add_argument('--dis_min',dest="dis_min",metavar='',default=None,type=int,\
                        help="Minimal distance between nucleosome dyads.")
    parser.add_argument('--dis_max',dest="dis_max",metavar='',default=None,type=int,\
                        help="Maximal distance between nucleosome dyads.")
    parser.add_argument('--dis_step',dest="dis_step",metavar='',default=10,type=int,\
                        help="Bin size or step used to calculate the distribution of nucleosome dyads distances.")
    
    parser.add_argument('--occ_min',dest="occ_min",metavar='',default=None,type=int,\
                        help="Minimal occupancy value of nucleosome peaks.")
    parser.add_argument('--occ_max',dest="occ_max",metavar='',default=None,type=int,\
                        help="Maximal occupancy value of nucleosome peaks.")
    parser.add_argument('--occ_step',dest="occ_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate occupancy values distribution.")
    
    parser.add_argument('--fuz_min',dest="fuz_min",metavar='',default=None,type=int,\
                        help="Minimal fuzziness score of nucleosome peaks.")
    parser.add_argument('--fuz_max',dest="fuz_max",metavar='',default=None,type=int,\
                        help="Maximal fuzziness score of nucleosome peaks.")
    parser.add_argument('--fuz_step',dest="fuz_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate fuzziness scores distribution.")
    
    parser.add_argument('--plot_row',dest="plot_row",metavar='',default=2,type=int,\
                        help="Number of rows to plot on each page.")
    parser.add_argument('--plot_column',dest="plot_column",metavar='',default=2,type=int,\
                        help="Number of columns to polt on each page.")
    parser.add_argument('--plot_colors', dest="plot_colors",metavar='',default='black,gray,red,blue,orange,purple,skyblue,cyan,green,blue4,darkgoldenrod',\
                        help="The colors to be used in the plot.")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print "\nfor more help, please try: python dantools peaksAnalysis -h\n"
            return 0
    else:
        print "\nfor help, please try: python dantools peaksAnalysis -h\n"
        return 0 
    
    print '\ncommand:\npython'," ".join(sys.argv) # print the command line, this let the user to keep a log and remember what parameters they specified

    occ,fuz={},{}
    peakfiles=args.peakfile_paths.split(',')
    peakaliases=args.peakfile_aliases.split(',')
    if len(peakfiles)!=len(peakaliases):
        print '\nError: peak files and aliases counts are not the same!'
        print len(peakfiles),'peak files:',peakfiles
        print len(peakaliases),'alias:',peakaliases
        print 'please note that aliases for different files would be separated by the comma \',\'\n'
        return 0
    for i in range(len(peakfiles)):
        peakfile,filealias=peakfiles[i],peakaliases[i]
        aliaspair=filealias.split(':')
        title=open(peakfile).readline().split()
        if len(aliaspair)>1:
            if 'smt_pos' in title:
                print '\nWrong:\nMore than one alias are provided for',peakfile,':',aliaspair
                print 'Please note that each peak file pooled/*peaks.xls (generated by DANPOS) would have only one alias.\n'
                return 0
            elif not ( ('fuzziness_diff_FDR' in title) and ('0-log10fuzziness_diff_pval' in title) ):
                print '\nWrong:\nThe input file',peakfile,'is not generated by DANPOS-2.1.1 or later version.'
                print 'If your input files are generated by DANPOS-2.1.0 or earlier version, please use the files *peaks.xls under the directory \"pooled/\".\n'
                return 0
            occ[aliaspair[1]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='control_smt_loca',val_col_name='control_smt_val',direction_by=[],top_value=None,bottom_value=None)
            occ[aliaspair[0]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='treat_smt_loca',val_col_name='treat_smt_val',direction_by=[],top_value=None,bottom_value=None)
            fuz[aliaspair[1]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='control_smt_loca',val_col_name='control_fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
            fuz[aliaspair[0]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='treat_smt_loca',val_col_name='treat_fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
        else:
            if not 'fuzziness_score' in title:
                print '\nWrong:\nThe input file',peakfile,"doesn't contain the column 'fuzziness_score'."
                print "If your input file is generated by DANPOS-2.1.1, please use the file *allpeaks.xls and specify a pair of alias for the control sample and treat sample.\n"
                return 0
            occ[aliaspair[0]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='smt_pos',val_col_name='smt_value',direction_by=[],top_value=None,bottom_value=None)
            fuz[aliaspair[0]]=retrieve_peaks_by_value(in_file=peakfile,out_file=None,cr_col_name='chr',pos_col_name='smt_pos',val_col_name='fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
    print ''
    
    if args.occ_max==None or args.occ_min==None:
        minmax=peakDicMinMax(occ,lowPercent=1,highPercent=99)
        if args.occ_min==None:args.occ_min=minmax[0]
        if args.occ_max==None:args.occ_max=minmax[1]
    if args.occ_step==None:args.occ_step=(args.occ_max-args.occ_min)/100.0
        
    if args.fuz_max==None or args.fuz_min==None:
        minmax=peakDicMinMax(fuz,lowPercent=1,highPercent=99)
        if args.fuz_min==None:args.fuz_min=minmax[0]
        if args.fuz_max==None:args.fuz_max=minmax[1]
    if args.fuz_step==None:args.fuz_step=(args.fuz_max-args.fuz_min)/100.0
    
    if args.dis_min==None:args.dis_min=100
    if args.dis_max==None:args.dis_max=250
    if args.dis_step==None:args.dis_step=1
    
    disdic=batchPeakDistanceDistribution(occ,outname=args.name+'.peak_distance',min=args.dis_min,max=args.dis_max,step=args.dis_step)
    rcode=plot(dic=disdic,outname='',main='peak_distance_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.dis_min,xmax=args.dis_max,ymin=None,ymax=None,xlab='Distance',ylab='Percent',colors=args.plot_colors.split(','))
    occdic=batchPeakValDistribution(occ,outname=args.name+'.peak_occupancy',min=args.occ_min,max=args.occ_max,step=args.occ_step)
    rcode+=plot(dic=occdic,outname='',main='peak_occupancy_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.occ_min,xmax=args.occ_max,ymin=None,ymax=None,xlab='Occupancy',ylab='Percent',colors=args.plot_colors.split(','))
    fuzdic=batchPeakValDistribution(fuz,outname=args.name+'.peak_fuzziness',min=args.fuz_min,max=args.fuz_max,step=args.fuz_step)
    rcode+=plot(dic=fuzdic,outname='',main='peak_fuzziness_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.fuz_min,xmax=args.fuz_max,ymin=None,ymax=None,xlab='Fuzziness',ylab='Percent',colors=args.plot_colors.split(','))

    rcode='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'+rcode
    rcode='pdf("'+args.name+'.pdf")\n'+rcode
    rcode+='dev.off()\n'
    fo=open(args.name+'.R','w')
    fo.write(rcode)
    fo.close()
    r(rcode)
    print ''
    

def peakDicMinMax(dic,lowPercent=0,highPercent=100):
    outmin,outmax,values=0,0,[]
    for name in dic:
        for chr in dic[name]:
            values+=dic[name][chr].values()
            #minv,maxv=min(values),max(values)
            #if outmin>minv:outmin=minv
            #if outmax<maxv:outmax=maxv
    values.sort()
    lth=len(values)
    outmin,outmax=values[int(lth*lowPercent/100.0)],values[int(lth*highPercent/100.0)]
    return [outmin,outmax]

def runPeakSelector():
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the function peakSelector, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\npython dantools.py peakSelector  <peakfile_path> [optional arguments]\n\nfor more help, please try: python dantools peakSelector -h\n"
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython dantools.py peakSelector <peakfile_path>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="This tells dantools to run peakSelector")

    parser.add_argument('peakfile_path',default=None,\
                        help="A path to the peak file, should be in the default output format of DANPOS (the file *.allpeaks.xls or a file *.peaks.xls under the directory 'pooled/'). \
                        See DANPOS documentation for more information.")

    parser.add_argument('--valueSelector', dest="valueSelector",metavar='',default=None,\
                        help="Selection requirements defined in the format: 'name1:50:1000,name2::1000,and', \
                        'name1' and 'name2' are column names in the peak file, \
                        e.g. 'treat_smt_val:50:1000' means selecting treat_smt_val value between 50 than 1000, \
                        'control_smt_val::1000' means selecting control_smt_val value less than 1000, \
                        the 'and' at the end means all selections are required, relace 'and' by 'or' \
                        if only and at least one or selections is required.")
    parser.add_argument('--genicSelector', dest="genicSelector",metavar='',default=None,\
                        help="One or several names selected among transcription start site (TSS), \
                        transcription terminal site (TTS), coding start site (CSS), coding terminal site (CTS), exon start site (ESS), \
                        or exon terminal site (ETS).  Each name\
                        comes with a pair of head and end locations in the format: 'site1:head:end,site2:head:end,and', \
                        e.g. 'TSS:-350:50' means selecting peaks that are in the TSS flanking region \
                        from 350bp upstream to 50bp downstream. The 'and' at the end means all selections are required, relace 'and' by 'or'\
                        if only and at least one selections is required.")

    parser.add_argument('--peak_out', dest="peak_out",metavar='',default='peakSelector.out.xls',\
                        help="An output file name for saving the selected peaks.")

    parser.add_argument('--gene_file', dest="gene_file",metavar='',default=None,\
                        help="A reference gene set, required when need to select nucleosomes by genicSelector. \
                        Gene set file must contain at least the following columns ordered as : name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, \
                        we suggest to download gene file from the UCSC tables at http://genome.ucsc.edu/cgi-bin/hgTables?command=start ")
    '''
    parser.add_argument('--gene_out', dest="gene_out",metavar='',default=None,\
                        help="A name of the file for saving the genes associated with selected peaks, as defined by genicSelector")
    '''
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0
    
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print "\nfor more help, please try: python dantools peakSelector -h\n"
            return 0
    else:
        print "\nfor help, please try: python dantools peakSelector -h\n"
        return 0 
    
    print '\ncommand:\npython'," ".join(sys.argv) # print the command line, this let the user to keep a log and remember what parameters they specified
    retr=open(args.peakfile_path).readlines()
    if args.valueSelector!=None:retr=peakSelectorByValue(peakLines=retr,selection=args.valueSelector)
    print max(len(retr)-1,0),'peaks selected by values'
    if args.genicSelector!=None:
        retr=peakSelectorByGeneStructure(peakLines=retr,selection=args.genicSelector,geneFile=args.gene_file,outGeneFile=None,chrbinsize=1)
        print max(len(retr)-1,0), 'peaks selected by genic structure'
    if len(retr)>1:
        fo=open(args.peak_out,'w')
        if args.peak_out[-3:]=='bed':
            n=0
            for line in retr[1:]:
                n+=1
                col=line.split()
                if col[1]=='-' and col[2]=='-':col[1],col[2]=str(int(col[3])-74),str(int(col[3])+74)
                elif col[1]=='-':col[1]=str(int(col[2])-74)
                elif col[2]=='-':col[2]=str(int(col[1])+74)
                if int(col[1])>int(col[2]):
                    temp=col[1]
                    col[1]=col[2]
                    col[2]=temp
                fo.write('\t'.join(col[:3]+[col[0]+':'+col[1]+'-'+col[2],'0','+'])+'\n')
        else:
            for line in retr:fo.write(line)
        fo.close()
        print ''
    else:
        print 'No peaks selected!\n'
    
def retrievePeakValuesAtRanks():
    '''
    Description:
        None
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\npython dantools.py retrievePeakValuesAtRanks <peak_file_path>  [optional arguments] \n\nfor more help, please try: python dantools retrievePeakValuesAtRanks -h\n"
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython dantools.py retrievePeakValuesAtRanks <peakfile_path>  [optional arguments] \n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="This tells dantools to run retrievePeakValuesAtRanks")

    parser.add_argument('peakfile_path',default=None,\
                        help="a path to the peak file, should be in the default output format of DANPOS (the file *.allpeaks.xls or a file *.peaks.xls under the directory 'pooled/'). \
                        See DANPOS documentation for more information.")

    parser.add_argument('--valueRanks', dest="valueRanks",metavar="",default=None,\
                        help="Specify the column names and ranks of values which are to be returned, \
                        e.g. control_smt_val:50:100:-200,treat_smt_val:500:-100 will return the \
                        control_smt_val value ranked at 50 and 100 by increasing order and 200 by decreasing order, as well as\
                        the treat_smt_val value ranked at 500 by increasing order and 100 by decreasing order.\
                        The default is to return value at each quarter for each column name.")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0
    
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print "\nfor more help, please try: python dantools retrievePeakValuesAtRanks -h\n"
            return 0
    else:
        print "\nfor help, please try: python dantools retrievePeakValuesAtRanks -h\n"
        return 0 
    print ''
    peakLines=open(args.peakfile_path).readlines()
    tcol=peakLines[0].split()
    if args.valueRanks==None:
        #print 'Please specify ranks for one or several of the follow column names:'
        #print '\n'.join(tcol[1:])
        #print ''
        #return
        tlen=len(peakLines)-1
        args.valueRanks=''
        for i in range(4,len(tcol)):args.valueRanks+=tcol[i]+':1:'+str(tlen/4)+':'+str(tlen/2)+':'+str(tlen*3/4)+':'+str(tlen)+','
        args.valueRanks=args.valueRanks[:-1]
    sels=args.valueRanks.split(',')
    for tsel in sels:
        sel=tsel.split(':')
        colid=0
        if not sel[0] in tcol[1:]:
            print "error:", sel[0],'is not a column name in the peak file'
            return 0
        else:
            for i in range(len(tcol)):
                if tcol[i]==sel[0]:colid=i
            if colid==0:
                print 'error: can not rank by', tcol[colid]
                return 0
        
        print sel[0]
        values=[]
        for line in peakLines[1:]:
            col=line.split()
            try:values.append(float(col[colid]))
            except:print col
        values.sort()
        for r in sel[1:]:
            r=int(r)
            if r>0:print '\trank',r,':',values[r-1]
            else:print '\trank',r,':',values[r]
    print ''
    
if __name__ == "__main__":
    if len(sys.argv)>1:
        if sys.argv[1]=='wigAnalysis':runWiggleAnalysis()
        elif sys.argv[1]=='peaksAnalysis':runPeaksAnalysis()
        elif sys.argv[1]=='peakSelector':runPeakSelector()
        elif sys.argv[1]=='retrievePeakValuesAtRanks':retrievePeakValuesAtRanks()
        else:printHelp()
    else:printHelp()

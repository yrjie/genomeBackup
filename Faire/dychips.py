#!/usr/bin/env python
#Version 1.1.0

import os,sys
from time import time
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow dychips to print each message on screen immediately.

def rundychips():
    '''
    Description:
        this function provide an entrance to the package dychips.
        It parse input parameters, print some help messages if required, else it call and pass all parameters values to the main function dychips().
    
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython dychips.py <path> [optional arguments]\n\nfor more help, please try: python dychips.py -h\n"
        return 0
    
    # parse all input parameters
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\n\npython dychips.py <path> [optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    parser.add_argument('path',default=None,\
                        help="Pairs of paths to ChIPseq data sets, \
                        e.g. file1.bed:dir2/,dir3/, \
                        the two paths in each pair must be seperated by ':', \
                        different pairs must be seperated by ',', \
                        each path could be directing to a file or a directory containing multiple files, \
                        with each file contains the data for one replicate, and files under one directory \
                        represent multiple replicates for the same group.")
    parser.add_argument('-b','--bg',dest='bg',metavar='',default=None,\
                        help="pairs of paths, each pair secify a genomic background data set for a ChIPseq data set, \
                        e.g. file1.bed:bgdir1,dir2/:None,dir3/:bg3.bed\
                        put a word 'None' when a ChIPseq data set has no background data set \
                        ")
    parser.add_argument('-c','--count',dest='count',metavar='',default=None,\
                        help="expected count of reads per replicate for each ChIPseq data set, e.g. \
                        file1.bed:10000000,dir2/:20000000,dir3/:15000000")
    parser.add_argument('-o','--out', dest="name",metavar='',default='result',\
                        help="a name for the output directory")
    parser.add_argument('-p', '--poccu',dest="poccu",metavar='',default=1e-5,type=float,help="P value cutoff for occupancy peak calling")
    parser.add_argument('-q', '--height',dest="height",metavar='',default=0,type=int,\
                        help="occupancy/intensity cutoff for peak calling")
    parser.add_argument('-t', '--testcut',dest="testcut",metavar='',default='1e-5',\
                        help="P value cutoff for calling differential peaks, or foldchange cutoff when -s is set to F")
    parser.add_argument('-n', '--nor',dest="nor",metavar='',default='F',\
                        help="data normalization method, could be 'F','Q','S' or 'N', \
                        representing quantile normalization, fold normalization, normalize by sampling, or no normalization")
    parser.add_argument('-w', '--width',dest="width",metavar='',default=40,type=int,\
                        help="minimal width of peak")
    parser.add_argument('-d', '--distance',dest="distance",metavar='',default=500,type=int,\
                        help="minimal distance between peaks, peaks closer than d will be merged as one peak")
    #parser.add_argument('-e', '--edge',dest="edge",metavar='',default=0,type=int,\
    #                    help="set to 1 if need to detect edges for occupancy/intensity peaks,else set to 0")
    parser.add_argument('-a', "--span", dest="span",metavar='',default=10,type=int,\
                        help="the span or step size  in the generated wiggle data")
    parser.add_argument('-k', '--keep',dest="keep",metavar='',default=0,type=int,\
                        help="save middle stage files? set to 0 if don't save, otherwise set to 1")
    parser.add_argument('-z', '--smooth_width', dest="smooth_width",metavar='',default=20,type=int,\
                        help="the smooth width before peak calling, set to 0 if need not to smooth")
    parser.add_argument('-x', "--pcfer", dest="pcfer",metavar='',default=0,type=int,\
                        help="set to 1 if want to do peak calling for each replicate,else set to 0")
    parser.add_argument('-l', '--lmd',dest="lmd",metavar='',default=300,\
                        help="lambda width for smoothing background data before background subtraction")
    parser.add_argument('-s', '--statis',dest="statis",metavar='',default='P',\
                        help="the statistics method for differential test, could be 'C','P', 'F', or 'S', \
                        representing Chi-square test, Possion test, Fold change, or direct subtraction")
    #parser.add_argument('-g', '--gapfill',dest="gapfill",metavar='',default=0,type=int,\
    #                    help="do gap filling? for peak calling, fill gap between peaks \
    #                    if the gap size is similar to peak size, set to 0 if don't fill, otherwise set to 1")
    parser.add_argument("--paired", dest="paired",metavar='',default=0,type=int,\
                        help="# ignore this when '-i' is 'wig', set to 1 if the data is paired-end reads")
    parser.add_argument("--frsz", dest="fs",metavar='',default=None,type=int,\
                        help="# ignore this when '-i' is 'wig',average size of \
                        the input DNA fragments in the seuqnecing experiment")
    parser.add_argument("--mifrsz", dest="mifrsz",metavar='',default=50,type=int,\
                        help="# ignore this when '-i' is 'wig',estimated minimal \
                        size of the DNA fragments, effective only when '--frsz' is not provided")
    parser.add_argument("--mafrsz", dest="mafrsz",metavar='',default=250,type=int,\
                        help="# ignore this when '-i' is 'wig',estimated maximal size of the \
                        DNA fragments, effective only when '--frsz' is not provided")
    parser.add_argument("--extend", dest="extend",metavar='',default=80,type=int,\
                        help="# ignore this when '-i' is 'wig',expected average size of DNA fragments, \
                        the real average size will be extend to this size when reads data is transform into occupancy data")
    parser.add_argument("--clonalcut", dest="clonalcut",metavar='',default='1',\
                        help="the cutoff for adjusting clonal signal, \
                        set as a P value larger than 0 and smaller than 1,\
                        set as 0 if don't need to adjust clonal signal,\
                        or set as 1 to allow automative detection.")
    
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print "\ndychips version 1.1.0\n"
        parser.print_help()
        print "\n"
        return 0
    elif len(sys.argv)>=2: # at least one parameter need to be specified: a path to input file or files
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print "\nfor more help, please try: python dychips.py -h\n"
            return 0

    print "\ndychips version 1.1.0\n"
    print 'command:\npython'," ".join(sys.argv) # print the command line, this let the user to keep a log and remember what parameters they specified
    from functions import dychips  
    print '\n',args # print all parameter values, this provide a eacy way for the user to double check the parameter values used by dychips
    
    from math import log10
    if args.statis=='C' or args.statis=='P':
        temptestcut=args.testcut.split('-')
        if len(temptestcut)==2:testcut=float(temptestcut[1])-log10(float(temptestcut[0][:-1]))
        else:testcut=0-log10(float(temptestcut[0]))
    elif args.statis=='F':testcut=log(float(args.statis))/log(2)
    
    dychips(opath=args.name,tpath=args.path,tbg=args.bg,\
           amount=args.count,edge=1,nor=args.nor,\
           test=args.statis,save=args.keep,both=True,\
           logp=testcut,width=args.width,distance=args.distance,\
           fill_gap=0,fill_value=1,pheight=args.poccu,height=args.height,\
           lmd=args.lmd,fs=args.fs,cut=args.clonalcut,\
           wgfmt='fixed',step=args.span,extend=args.extend,\
           mifrsz=args.mifrsz,mafrsz=args.mafrsz,pcfer=args.pcfer,\
           smooth_width=args.smooth_width,paired=args.paired)                  

if __name__ == "__main__":
    rundychips()

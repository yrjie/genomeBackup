import sys,os
from optparse import OptionParser

#parser=OptionParser(usage='''%prog''')
parser=OptionParser()
parser.add_option('-a','--asm',help='Assembly name', dest='asm')
#print asm
parser.add_option('-v',action='store_true',default=False, dest='verbose')
(opts, args)=parser.parse_args()
print opts, args
if opts.asm==None:
    parser.print_help()
    exit(1)

import sys,os
from subprocess import *

if len(sys.argv)<2:
    print 'Usage: python getInterNum.py windowSize'
    exit(1)

file1='CMH009_peaks.bed'
file2='wgEncodeLicrHistoneHeartH3k4me3MAdult8wksC57bl6StdPk.broadPeak'
file3='wgEncodeLicrHistoneHeartH3k04me3UE14halfC57bl6StdPk.broadPeak'

cmd='windowBed -a '+file1+' -b '+file2+' -w '+sys.argv[1]+' -u'
p=Popen(cmd,stdout=PIPE)
n12=len(p.stdout.readlines())
cmd='windowBed -a '+file1+' -b '+file3+' -w '+sys.argv[1]+' -u'
p=Popen(cmd,stdout=PIPE)
n13=len(p.stdout.readlines())
cmd='windowBed -a '+file2+' -b '+file3+' -w '+sys.argv[1]+' -u'
p=Popen(cmd,stdout=PIPE)
n23=len(p.stdout.readlines())

cmd='R --no-save --slave --args '+str(n12)+' '+str(n13)+' '+str(n23)+' <plotVenn.r'
os.system(cmd)

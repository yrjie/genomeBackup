#!/home/ruijie/bin/python

import os,sys
import numpy

if len(sys.argv)<2:
    print 'Usage: python plotCorr.py infile'
    exit(1)

fi=open(sys.argv[1])
dat=[]
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=[float(x) for x in line.split()]
    dat.append(temp)
fi.close()
corr=numpy.corrcoef(dat)
num=corr.shape[0]
mcorr=str((corr.sum()-num)/num/(num-1))
print 'Mean corr: '+mcorr
numpy.savetxt('corr.tab',corr,fmt='%.4f',delimiter='\t')
plotC='matrix2png -data corr.tab -mincolor blue -maxcolor green -size 1:1 -s >corr.png -title "MeanCorr: '+mcorr+'"'
os.system(plotC)

import sys,os
import statsmodels.api as sm

if len(sys.argv)<3:
    print 'Usage: dat1 dat2'
    exit(1)

eps=10
dat=[[] for i in range(2)]
ecdf=[[] for i in range(2)]
for i in range(2):
    fi=open(sys.argv[i+1])
    for line in fi:
    	line=line.strip()
	if len(line)<1:
	    continue
	dat[i].append(float(line))
    ecdf[i]=sm.distributions.ECDF(dat[i])
    fi.close()

mi=min([min(dat[0]),min(dat[1])])
ma=max([max(dat[0]),max(dat[1])])
for i in range(11):
    now=mi+(ma-mi)*i/10.0
    #pr0=ecdf[0](now)-ecdf[0](now-eps)
    #pr1=ecdf[1](now)-ecdf[1](now-eps)
    pr0=ecdf[0](now)
    pr1=ecdf[1](now)
    print '\t'.join([str(x) for x in [now,pr0, pr1, pr0-pr1]])

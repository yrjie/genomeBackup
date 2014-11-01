import sys,os

# the 4th field should be the identifier

if len(sys.argv)<2:
    print 'Usage: overlap'
    exit(1)

fi=open(sys.argv[1])
now=''
one=150
binN=10
num=[1]*binN
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[3]==now:
    	leng=int(temp[-2])
	if leng<0:
	    leng=-leng
    	idx=leng/one
	if idx>=binN:
	    idx=binN-1
	num[idx]+=1
    else:
	if now!='':
	    print '\t'.join([str(1.0*x/sum(num)) for x in num])
    	now=temp[3]
	num=[1]*binN
print '\t'.join([str(1.0*x/sum(num)) for x in num])
fi.close()

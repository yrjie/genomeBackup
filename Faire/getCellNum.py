import sys,os

if len(sys.argv)<2:
    print 'Usage: combinedPeak'
    exit(1)

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    num=len(set(temp[3].split(';')))
    if num<2:
    	continue
    temp.append(str(num))
    print '\t'.join(temp)
fi.close()

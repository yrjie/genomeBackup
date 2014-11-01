import sys,os

if len(sys.argv)<2:
    print 'Usage: python chop2bw.py bam.lst'
    exit(1)

fi=open(sys.argv[1],'r')
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('.')
    print temp[0]
    os.system('sh runChop.sh '+temp[0])
    os.system('sh runToBW.sh '+temp[0])
fi.close()

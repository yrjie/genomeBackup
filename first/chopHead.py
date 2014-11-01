import sys,os

if len(sys.argv)<2:
    print 'Usage: python chopHead.py bedFile'
    exit(1)

fi=open(sys.argv[1],'r')
fo=open('out.bed','w')
for line in fi:
    temp=line.strip().split('\t')
    if len(temp)<6:
    	continue
    if temp[5]=='+':
    	temp[2]=str(int(temp[1])+2)
    else:
    	temp[1]=str(int(temp[2])-2)
    for j in range(0,6):
    	fo.write(temp[j])
	if j<5:
	    fo.write('\t')
	else:
	    fo.write('\n')
fi.close()
fo.close()


import sys,os

if len(sys.argv)<2:
    print 'Usage: file(2 or 3 lines) [name1] [name2] [name3]'
    exit(1)

nameA='A'
nameB='B'
nameC='C'
if len(sys.argv)>=4:
    nameA=sys.argv[2]
    nameB=sys.argv[3]
if len(sys.argv)>=5:
    nameC=sys.argv[4]
fi=open(sys.argv[1])
tempA=fi.readline().split('\t')
tempB=fi.readline().split('\t')
try:
    line=fi.readline().strip()
    if len(line)<1:
    	tempC=[]
    else:
    	tempC=line.split('\t')
except:
    tempC=[]
fi.close()
setA=set(tempA)
setB=set(tempB)
setC=set(tempC)
if len(setC)==0:
    num1=str(len(setA))
    num2=str(len(setA&setB))
    num3=str(len(setB))
    print '\t'.join([num1, num2, num3])
    cmd=' '.join(['R --no-save --slave --args',num1,num2,num3,nameA,nameB,'<~/bin/plotVenn.r'])
else:
    n1=str(len(setA))
    n2=str(len(setB))
    n3=str(len(setC))
    n12=str(len(setA&setB))
    n13=str(len(setA&setC))
    n23=str(len(setB&setC))
    n123=str(len(setA&setB&setC))
    print '\t'.join([n1,n2,n3,n12,n23,n123])
    cmd=' '.join(['R --no-save --slave --args ', n1, n2, n3, n12, n13, n23, n123, nameA, nameB, nameC, '<~/bin/plotVenn.r'])
os.system(cmd)

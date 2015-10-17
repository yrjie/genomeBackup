import sys,os

if len(sys.argv)<2:
    print 'Usage: peakTag [norm]'
    print 'output will have a key at the head'
    exit(1)

file=sys.argv[1]

binN=40
left=-500
right=500
binS=(right-left)/binN
allRatio=[]

norm=0
if len(sys.argv)>2:
    norm=int(sys.argv[2])

def getRatio(infile):
    nowId=''
    tagOne=[0]*(binN+1)
    fi=open(infile)
    for line in fi:
        line=line.strip()
        if len(line)<1:
            continue
        temp=line.split('\t')
        fragLen=int(temp[-2])
        id='\t'.join([temp[0],temp[1],temp[2]])
        if nowId=='':
            nowId=id
        if id != nowId:
            if tagOne[-1]>0:
                for i in range(binN):
                    tagOne[i]/=1.0*tagOne[-1]
                for i in range(binN/2):
                    tagOne[i]-=tagOne[binN-1-i]
                #ma=max(tagOne[0:binN])
		#for i in range(binN):
		#    tagOne[i]/=ma
#             print '\t'.join([nowId] + [str(x) for x in tagOne])
#             print '\t'.join([str(x) for x in tagOne])
            allRatio.append(tagOne)
            tagOne=[0]*len(tagOne)
            nowId=id
        if fragLen<left or fragLen>=right:
            continue
        tagOne[(fragLen-left)/binS]+=1
        tagOne[-1]+=1
    fi.close()
    for i in range(binN):
        tagOne[i]/=1.0*tagOne[-1]
    allRatio.append(tagOne)
#     print '\t'.join([nowId] + [str(x) for x in tagOne])
#     print '\t'.join([str(x) for x in tagOne])

def printOut(norm):
    meanR=[0]*(binN)
    for x in allRatio:
        for i,y in enumerate(meanR):
            meanR[i]+=x[i]*x[-1]
    numTag=sum(meanR)
    for i,y in enumerate(meanR):
        meanR[i]/=numTag
    meanR.append(1)
    if norm:
        for id,x in enumerate(allRatio):
            #print '\t'.join([str(id)]+[str(y/meanR[i]) for i,y in enumerate(x)])
            print '\t'.join([str(y/meanR[i]) for i,y in enumerate(x)])
    else:
        for id,x in enumerate(allRatio):
            #print '\t'.join([str(id)]+[str(y) for i,y in enumerate(x)])
            print '\t'.join([str(y) for i,y in enumerate(x)])

getRatio(file)
printOut(norm)

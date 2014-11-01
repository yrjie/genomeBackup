import sys,os

if len(sys.argv)<3:
    print 'Usage: python getPatternCount.py dnaseHead chipseqHead'
    exit(1)

fo=open('pattern/'+sys.argv[1]+'.xls','w')
fo.write('down/-\tdown/+\tup/-\tup/+\n')
cmd='bigWigSummaryBatch 2bp/'+sys.argv[1]+'.minus.bw data/'+sys.argv[2]+'_minus.narrowPeak 1 -type=max >temp0.lst'
os.system(cmd)
cmd='bigWigSummaryBatch 2bp/'+sys.argv[1]+'.plus.bw data/'+sys.argv[2]+'_minus.narrowPeak 1 -type=max >temp1.lst'
os.system(cmd)
cmd='bigWigSummaryBatch 2bp/'+sys.argv[1]+'.minus.bw data/'+sys.argv[2]+'_plus.narrowPeak 1 -type=max >temp2.lst'
os.system(cmd)
cmd='bigWigSummaryBatch 2bp/'+sys.argv[1]+'.plus.bw data/'+sys.argv[2]+'_plus.narrowPeak 1 -type=max >temp3.lst'
os.system(cmd)
data=[]
for i in range(4):
    fi=open('temp'+str(i)+'.lst','r')
    data.append(fi.readlines())
    fi.close()
for i in range(len(data[0])):
    for j in range(4):
    	data[j][i]=data[j][i].strip()
	if len(data[j][i])<1:
	    continue
	fo.write(data[j][i])
	if j==3:
	    fo.write('\n')
	else:
	    fo.write('\t')
fo.close()

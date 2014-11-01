import sys,os

if len(sys.argv)<2:
    print 'Usage: python addName.py tfbs.lst [span]'
    exit(1)

span=100
if len(sys.argv)>=3:
    span=int(sys.argv[2])

inDir='/home/browser/BASIC/import-data/narrowPeak/'
outDir='./data/'
fi=open(sys.argv[1],'r')

for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    t1=line.split('\t')
    file=inDir+t1[0]+'.narrowPeak'
    cell=line.split('\t')[1].split(' ')[0]
    fiData=open(file,'r')
    outfile=outDir+os.path.basename(file)
    print outfile
    fo=open(outfile,'w')
    for row in fiData:
    	row=row.strip()
	if len(row)<1:
	    continue
	temp=row.split('\t')
	if int(temp[9])<0:
	    peak=(int(temp[1])+int(temp[2]))/2
	else:
	    peak=int(temp[1])+int(temp[9])
	temp[1]=str(peak-span)
	temp[2]=str(peak+span)
	temp[3]=cell
	fo.write('\t'.join(temp)+'\n')
    fo.close()
    fiData.close()
fi.close()


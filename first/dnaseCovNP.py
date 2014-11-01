import sys,os

if len(sys.argv)<3:
    print 'Usage: python dnaseCovNP.py DNase.narrowPeak chipFileList [outfile]\nIf outfile is not provided, <prefix>.cov will be used'
    exit(1)

chipPath='/home/browser/BASIC/import-data/narrowPeak/'
temp=sys.argv[1].split('/')
outfile=temp[len(temp)-1].split('.')[0]+'.cov'
if len(sys.argv)==4:
    outfile=sys.argv[3]
print outfile

fi=open(sys.argv[2],'r')
fo=open(outfile,'w')

for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    line=chipPath+line
    pNum=len(open(line,'r').readlines())
    cmd='bedtools intersect -a '+sys.argv[1]+' -b '+line+' >temp'
#    print cmd
    os.system(cmd)
    cov=len(open('temp','r').readlines())
    cmd='bedtools shuffle -i '+line+' -g ~/genome/human.hg19.genome >shuffled'
#    print cmd
    os.system(cmd)
    cmd='bedtools intersect -a '+sys.argv[1]+' -b shuffled >tempS'
#    print cmd
    os.system(cmd)
    FP=len(open('tempS','r').readlines())
    temp=line.split('/')
    line=temp[len(temp)-1].split('.')[0]
    fo.write(str(pNum)+'\t'+str(1.0*cov/pNum)+'\t'+str(1.0*FP/pNum)+'\t'+line+'\n')

fi.close()
fo.close()

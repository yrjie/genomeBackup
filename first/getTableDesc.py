import sys,os

if len(sys.argv)<2:
    print 'Usage: python getTableDesc.py type.lst'
    exit(1)

fi=open(sys.argv[1],'r')
cnt=0
for line in fi:
    temp=line.strip().split('\t')
    outfile='desc'+str(cnt)
    fo=open(outfile,'w')
    fo.write('\n'+temp[0]+'\n')
    fo.close()
    cmd='mysql -uroot hg19 -e \"describe '+temp[1]+'\" >>'+outfile
    os.system(cmd)
    cnt+=1
fi.close()
cmd='cat'
for i in range(0,cnt):
    cmd+=' desc'+str(i)
os.system(cmd)

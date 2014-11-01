import os,sys

if len(sys.argv)<2:
    print 'Usage: swap.bedpe'
    exit(1)

fi=open(sys.argv[1])
foS=open('sample.bed','w')
foC=open('control.bed','w')
cnt=0
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    beg=int(temp[1])
    end=int(temp[5])
    leng=end-beg
    if leng<100:
    	foC.write('\t'.join([temp[0],temp[1],temp[5], str(cnt), '0', '+'])+'\n')
    	foC.write('\t'.join([temp[0],temp[1],temp[5], str(cnt), '0', '-'])+'\n')
    elif leng>180 and leng<247:
    	foS.write('\t'.join([temp[0], str(beg+50), str(end-50), str(cnt), '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(beg+50), str(end-50), str(cnt), '0', '-'])+'\n')
    elif leng>315 and leng<473:
    	foS.write('\t'.join([temp[0], str(beg+50), str(beg+150), str(cnt)+'.1', '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(end-150), str(end-50), str(cnt)+'.2', '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(beg+50), str(beg+150), str(cnt)+'.1', '0', '-' ])+'\n')
    	foS.write('\t'.join([temp[0], str(end-150), str(end-50), str(cnt)+'.2', '0', '-'])+'\n')
    elif leng>558 and leng<615:
    	foS.write('\t'.join([temp[0], str(beg+50), str(beg+150), str(cnt)+'.1', '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(beg+250), str(beg+350), str(cnt)+'.2', '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(end-150), str(end-50), str(cnt)+'.3', '0', '+'])+'\n')
    	foS.write('\t'.join([temp[0], str(beg+50), str(beg+150), str(cnt)+'.1', '0', '-'])+'\n')
    	foS.write('\t'.join([temp[0], str(beg+250), str(beg+350), str(cnt)+'.2', '0', '-'])+'\n')
    	foS.write('\t'.join([temp[0], str(end-150), str(end-50), str(cnt)+'.3', '0', '-'])+'\n')
    cnt+=1
fi.close()
foC.close()
foS.close()

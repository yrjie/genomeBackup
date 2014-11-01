import os
import sys
import glob


dirs=os.listdir('.')
for dir in dirs:
	if os.path.isdir(dir) == False:
		continue
        files=glob.glob(dir+'/*')
        outf=open(dir+'.txt','w')
        print(dir)
        for fl in files:
                fsize=int(os.path.getsize(fl))
                lines=open(fl,'r').readlines()
                fsize=fsize-len(lines)-len(lines[0])
                outf.write(fl.replace('.fa','').replace(dir+'/','')+'\t'+str(fsize)+'\n')
        outf.close()


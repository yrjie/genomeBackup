#!/home/ruijie/bin/python

import sys,os

if len(sys.argv)<4:
    print 'python tab2arff.py positive negative name'
    exit(1)

print '@relation '+sys.argv[3]+'\n'

fi=open(sys.argv[1],'r')
cnt=1
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if cnt==1:
	for i in range(1,len(temp)+1):
	    print '@attribute a'+str(i)+' numeric'
	print '@attribute class {p, n}\n\n@data\n'
    print ','.join(temp)+',p'
    cnt+=1
fi.close()

fi=open(sys.argv[2],'r')
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    print ','.join(temp)+',n'
fi.close()

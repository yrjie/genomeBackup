import sys,os

if len(sys.argv)<3:
    print 'Usage: python filterFaln.py bedfile threshold'
    exit(1)

#class Read:
#    chr=''
#    start=''
#    end=''
#    def __init__(self,chr,start,end):
#    	self.chr=chr
#    	self.start=start
#    	self.end=end

def overlap(line1,line2):
    temp1=line1.split('\t')
    temp2=line2.split('\t')
    return int(temp2[1])<=int(temp1[2])

lines=open(sys.argv[1]).readlines()
offset=int(sys.argv[2])
head=0
nowChr='chr0'
for i in range(len(lines)):
    line=lines[i].strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if nowChr!=temp[0]:
    	head=i
    	nowChr=temp[0]
    if i-head<offset:
    	print line
    	continue
    if not overlap(lines[head], line):
    	print line
    	head+=1

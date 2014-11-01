import sys,os

if len(sys.argv)<2:
    print 'Usage: combineFile'
    exit(1)

lnk=['GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC',
	'GCTGTTAAGGACCGTACATCCGCCTTGGCCGTGACATTGACC',
	'GGTCAATGTCACCGTACATCCGCCTTGGCCGTCCTTAACAGC',
	'GGTCAATGTCACCGTACATCCGCCTTGGCCGTGACATTGACC',
	'GCTGTTAAGGACGGCCAAGGCGGATGTACGGTCCTTAACAGC',
	'GGTCAATGTCACGGCCAAGGCGGATGTACGGTCCTTAACAGC',
	'GCTGTTAAGGACGGCCAAGGCGGATGTACGGTGACATTGACC',
	'GGTCAATGTCACGGCCAAGGCGGATGTACGGTGACATTGACC'
	]

pr15={'pr1':'GACGCTCTTCCGATC',
	'pr2':'ACCGCTCTTCCGATC',
	'pr1Rev':'GATCGGAAGAGCGTC',
	'pr2Rev':'GATCGGAAGAGCGGT'
	}

#print len(lnk), len(pr15)
# head: start point of tag, tail: end point +1
fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1 or line[0]=='>':
    	continue
    x=''
    for l in lnk:
	if l in line:
	    x=l
	    break
    if x=='':
    	continue
    head1=0
    tail1=line.index(x)
    head2=tail1+len(x)
    tail2=len(line)
    if pr15['pr1'] in line:
    	head1=line.index(pr15['pr1'])+len(pr15['pr1'])
    if pr15['pr1Rev'] in line:
    	head1=line.index(pr15['pr1Rev'])+len(pr15['pr1Rev'])
    if pr15['pr2'] in line:
    	tail2=line.index(pr15['pr2'])
    if pr15['pr2Rev'] in line:
    	tail2=line.index(pr15['pr2Rev'])
    tag1=tail1-head1
    tag2=tail2-head2
    if tag1<0 or tag2<0:
    	continue
    print str(tag1)+'\t'+str(tag2)
fi.close()


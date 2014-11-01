import sys,os

if len(sys.argv)<2:
    print 'Usage: combined.fastq'
    exit(1)

def print_fq(id, tag, qua):
    print id
    print tag
    print '+'
    print qua

#lnk=['GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC',
#	'GCTGTTAAGGACCGTACATCCGCCTTGGCCGTGACATTGACC',
#	'GGTCAATGTCACCGTACATCCGCCTTGGCCGTCCTTAACAGC',
#	'GGTCAATGTCACCGTACATCCGCCTTGGCCGTGACATTGACC',
#	'GCTGTTAAGGACGGCCAAGGCGGATGTACGGTCCTTAACAGC',
#	'GGTCAATGTCACGGCCAAGGCGGATGTACGGTCCTTAACAGC',
#	'GCTGTTAAGGACGGCCAAGGCGGATGTACGGTGACATTGACC',
#	'GGTCAATGTCACGGCCAAGGCGGATGTACGGTGACATTGACC'
#	]

lnk=['GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC',
	'GGTCAATGTCACCGTACATCCGCCTTGGCCGTGACATTGACC',
	'GCTGTTAAGGACGGCCAAGGCGGATGTACGGTCCTTAACAGC',
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
while 1:
    l1=fi.readline().strip()
    if len(l1)<1:
    	break
    l2=fi.readline().strip()
    l3=fi.readline().strip()
    l4=fi.readline().strip()
    x=''
    for l in lnk:
	if l in l2:
	    x=l
	    break
    if x=='':
    	continue
    head1=0
    tail1=l2.index(x)
    head2=tail1+len(x)
    tail2=len(l2)
    if pr15['pr1'] in l2:
    	head1=l2.index(pr15['pr1'])+len(pr15['pr1'])
    if pr15['pr1Rev'] in l2:
    	head1=l2.index(pr15['pr1Rev'])+len(pr15['pr1Rev'])
    if pr15['pr2'] in l2:
    	tail2=l2.index(pr15['pr2'])
    if pr15['pr2Rev'] in l2:
    	tail2=l2.index(pr15['pr2Rev'])
    tag1=l2[head1:tail1]
    q1=l4[head1:tail1]
    tag2=l2[head2:tail2]
    q2=l4[head2:tail2]
    print_fq(l1+'_1', tag1, q1)
    print_fq(l1+'_2', tag2, q2)
#    tag1=tail1-head1
#    tag2=tail2-head2
#    if tag1<0 or tag2<0:
#    	continue
#    print str(tag1)+'\t'+str(tag2)
fi.close()


import sys,os
import random

for i in range(50,200):
    for j in range(1000):
    	num1=int(random.random()*i)
    	a='A'*num1+'GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC'+'A'*(i-num1)
	if random.random()<0.1:
	    a='321321321321GACGCTCTTCCGATC'+a
	elif random.random()<0.2:
	    a='ewewqewqGATCGGAAGAGCGTC'+a
	elif random.random()<0.3:
	    a+='ACCGCTCTTCCGATC8798789798798'
	elif random.random()<0.4:
	    a+='GATCGGAAGAGCGGT3213213weqw'
	print a

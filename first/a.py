import sys,os
import numpy as np
from numpy import *
#identity(3)
import scipy
import sys
#scipy.test(10)
#np.test()

print sys.path[0]
class A:
    a=1
    b='ab'
    def __init__(self,t1,t2):
    	self.a=t1
    	self.b=t2
    def __eq__(self,other):
    	return self.a==other.a

print sys.argv[1]
ans=''
lines=open('/home/yrjie/studio/temp','r').readlines()
fo=open('/home/yrjie/studio/temp','w')
print len(lines)
for line in lines:
    if 't' in line:
    	fo.write(line)
fo.close()
print ans


def main(str1, str2):
    print str1+str2
    a=[1,2,3,4]
    return a

def useless():
    print sys.argv[0]
    print sys.path[0]
    print np.identity(3)-np.identity(3)
    print "abc"
    word="123456"
    b=array([[1,2,3],[4,5,6]]);
    print b.shape
    c=array([[2],[3],[4]]);
    #d=[[6],[7]]
    print dot(b,c)
    print str(len(word))+word[3:]+word[:3]
    print 0
    a=int(raw_input("lu"));
    if a<0:
	print -1
    elif a==0:
	print 0
    else:
	print range(-2,-7,-3)


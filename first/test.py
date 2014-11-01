import sys,os

fi=open('in3','r')
line=fi.readline()
temp=line.split(' ')
for tf in temp:
    print tf.split('=')[1]
fi.close()

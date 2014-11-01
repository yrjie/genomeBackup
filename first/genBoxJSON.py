import sys,os

if len(sys.argv)<3:
    print 'Usage: python genBoxJSON.py infile title'
    exit(1)

fi=open(sys.argv[1],'r')
text='{\"chart\":{\"type\":\"boxplot\",\"renderTo\": \"container_boxplot\"},\n\"title\":{\"text\":\"'+sys.argv[2]+'\"},\n\"credits\":{\"enabled\":false},\n'
xAxis='\"xAxis\":{\"labels\":{\"rotation\":-90,\"align\":\"right\"},\"categories\":['
series='\"series\":[{\"name\":\"'+sys.argv[2]+'\",\"data\":['
i=0
for line in fi:
    line=line.strip()
    if len(line)<1:
    	continue
    temp=line.split('\t')
    if temp[5]=='0' or int(temp[0].split('.')[1])<=4:
    	continue
    if i>0:
    	xAxis+=','
    	series+=','
    xAxis+='\"'+temp[0]+'\"'
    series+='['+temp[1]+','+temp[2]+','+temp[3]+','+temp[4]+','+temp[5]+']'
    i+=1
fi.close()
print '# of TF: '+str(i)
series+=']}]\n'
xAxis+=']},\n'
text=text+xAxis+series+'}'
fo=open('/home/ruijie/public_html/data/box.json','w')
fo.write(text)
fo.close()

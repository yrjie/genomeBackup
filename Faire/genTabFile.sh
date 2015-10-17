if [ $# -lt 2 ]
then
    echo 'Usage: AtacPk tagfile'
    exit
fi

windowBed -a $1 -b $2 -w 0|awk 'BEGIN{pre="";
a[1]=100;a[2]=180;a[3]=247;a[4]=315;a[5]=473;a[6]=558;a[7]=615;a[8]=1000000;
for (i in a)
    cnt[i]=0
}
$9<0 {$9=-$9}
{now=$1$2$3;
    if (now!=pre){
    	if (pre!=""){
	    if (num==0) num=1
    	    for (i=1;i<=8;i++){
    		printf cnt[i]/num"\t";
    		cnt[i]=0;
    	    }
    	    printf num"\n"
    	    num=0
    	}
    	pre=now;
    }
    for (i=1;i<=8;i++){
	if ($9<a[i]){
	    cnt[i]++;
	    break;
	}
    }
    num++;
}'

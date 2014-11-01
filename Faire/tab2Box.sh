if [ $# -lt 2 ]
then
    echo 'Usage: tabfile name'
    exit
fi

awk -v name=$2 '
{
    for (i=1;i<=NF;i++){
	arr[NR,i]=$i;
	if(big <= NF)
	    big=NF;
    	}
}
	
END {
for(i=1;i<=big;i++){
    printf name"_"i;
    for(j=1;j<=NR;j++){
	printf("\t%s",arr[j,i]);
    }
    printf("\n");
}
}' $1 >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

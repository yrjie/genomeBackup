if [ $# -lt 3 ]
then
    echo "Usage: sh getRoc.sh peakfile cell suffix"
    exit
fi

file=$1
cell=$2
suffix=$3
rm output/$cell$suffix.tab
num=`wc -l $file |cut -d' ' -f1`
one=`echo $num/10 |bc`

#max1=10000 # for ATAC under DFilter
#max1=20000 # for ATAC
#max1=50000 # for ATAC_low
max1=100000 # for Faire
#max1=200000 # for Faire

if [ $max1 -lt $one ]
then
    one=$max1
fi
sortFile=sorted$RANDOM
# be careful about the column of score
sort -k4 -gr $file >$sortFile
#cp $file $sortFile
echo $cell$suffix.tab
for ((i=1;i<=10;i++)) do
    now=`echo $one*$i|bc`
    head -n $now $sortFile >tempRoc$suffix
    sh checkFpCov.sh $cell tempRoc$suffix $suffix
    awk 'BEGIN {FS=OFS="\t";fp=0;cov=0} {fp+=$5;cov+=$4} END {print fp/NR,cov/NR}' output/$cell$suffix'_all'.xls >>output/$cell$suffix.tab
done
rm output/$cell$suffix'_all'.xls
rm $sortFile

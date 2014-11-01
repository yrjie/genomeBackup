if [ $# -lt 3 ]
then
    echo "Usage: peakfile cell suffix [hmmPrefix]"
    exit
fi

file=$1
cell=$2
suffix=$3"Pk"
rm output/$cell$suffix.tab

prefix="pile"

if [ $# -gt 3 ]
then
    prefix=$4
fi

posPk=$prefix"hmmPosPk.bed"
negPk=$prefix"hmmNegPk.bed"

num=`wc -l $file |cut -d' ' -f1`
one=`echo $num/10 |bc`
posN=`wc -l $posPk|awk '{print $1}'`
negN=`wc -l $negPk|awk '{print $1}'`

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
    covPosN=`windowBed -a $posPk -b tempRoc$suffix -u -w 0|wc -l`
    covNegN=`windowBed -a $negPk -b tempRoc$suffix -u -w 0|wc -l`
#    echo $now $covPosN $covNegN
    TPR=`echo $covPosN/$posN|bc -l`
    FPR=`echo $covNegN/$negN|bc -l`
    echo -e $FPR"\t"$TPR >>output/$cell$suffix.tab
done

rm tempRoc$suffix $sortFile

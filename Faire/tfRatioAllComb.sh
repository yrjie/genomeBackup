if [ $# -lt 3 ]
then
    echo 'Usage: peakFile prefix tfFile'
    echo 'The 4th column is the signal'
    echo 'sort by the ratio of covering TF'
    exit
fi

prefix=$2
tfFile=$3
sig=(a b c d e f g h i j)
abcType=(ab abc ac aUnique bc bUnique cUnique)
tmpFile=tmp$RANDOM

sort -k4gr $1 >$tmpFile
num=`wc -l $1|awk '{print $1}'`
one=`echo $num/10|bc`
split $tmpFile -l $one split

for i in ${sig[*]}
do
    for j in ${abcType[*]}
    do
	afile=splita$i
	bfile=pileall3_$j.bed
	windowBed -a $afile -b $bfile -u -w 0 >$tmpFile
	num=`wc -l $tmpFile|awk '{print $1}'`
	echo $i $j $num
	numTF=`windowBed -a $tmpFile -b $tfFile -u -w 0|wc -l`
	#numTF=`windowBed -b $tmpFile -a $tfFile -u -w 0|wc -l`
	score=`echo $numTF/$num|bc -l`
	awk -v sc=$score 'BEGIN{OFS="\t"}{$4=sc;print $0}' $tmpFile >$prefix"_"$i"_"$j.bed
    done
done

awk -v sc=0 'BEGIN{OFS="\t"}{$4=sc;print $0}' splitak >$prefix"_"k.bed

cat $prefix* |sort -k4gr >$prefix.bed
rm $prefix"_"* $tmpFile splita*

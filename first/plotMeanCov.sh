if [ $# -lt 3 ]
then
    echo 'Usage: bwfile bedfile bin'
    exit
fi

bigWigSummaryBatch $1 $2 $3 >temp.dat
name=`awk '{if (NR==1) print "+-"(int(($3-$2)/2))}' $2`
num=`wc -l $2 |cut -d' ' -f1`
R --no-save --slave --args temp.dat $num:$name <~/bin/plotMeanCov.R

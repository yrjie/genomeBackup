if [ $# -lt 6 ]
then
    echo 'Usage: bed1 bed2 bwfileP bwfileM name1 name2'
    exit
fi

tmpfile="temp.bed"
win=200
awk -v win=$win '{print $1,$2-win,$3+win}' $1 >$tmpfile
bigWigSummaryBatch $3 $tmpfile 4 >tempp1
bigWigSummaryBatch $4 $tmpfile 4 >tempm1
#paste tempp1 tempm1 |awk -v name=$5"_L" 'BEGIN {printf name}{a=$1+0.1;b=$5+0.1;if (0&&a==b) next; printf "\t"log(a/b)}END{printf "\n"}' >temp.dat
#paste tempp1 tempm1 |awk -v name=$5"_R" 'BEGIN {printf name}{a=$4+0.1;b=$8+0.1;if (0&&a==b) next; printf "\t"log(a/b)}END{printf "\n"}' >>temp.dat
paste tempp1 tempm1 |awk '{a=$1+0.1;b=$5+0.1;print log(a/b)}' >temp1
paste tempp1 tempm1 |awk '{a=$4+0.1;b=$8+0.1;print log(a/b)}' >temp2
paste temp1 temp2 >pos.tab

awk -v win=$win '{print $1,$2-win,$3+win}' $2 >$tmpfile
bigWigSummaryBatch $3 $tmpfile 4 >tempp1
bigWigSummaryBatch $4 $tmpfile 4 >tempm1
#paste tempp1 tempm1 |awk -v name=$6"_L" 'BEGIN {printf name}{a=$1+0.1;b=$5+0.1;if (0&&a==b) next; printf "\t"log(a/b)}END{printf "\n"}' >>temp.dat
#paste tempp1 tempm1 |awk -v name=$6"_R" 'BEGIN {printf name}{a=$4+0.1;b=$8+0.1;if (0&&a==b) next; printf "\t"log(a/b)}END{printf "\n"}' >>temp.dat
paste tempp1 tempm1 |awk '{a=$1+0.1;b=$5+0.1;print log(a/b)}' >temp1
paste tempp1 tempm1 |awk '{a=$4+0.1;b=$8+0.1;print log(a/b)}' >temp2
paste temp1 temp2 >neg.tab
#R --no-save --slave --args temp.dat <~/bin/plotViolin.R

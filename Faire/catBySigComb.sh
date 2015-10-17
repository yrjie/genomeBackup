if [ $# -lt 1 ]
then
    echo 'Usage: peakFile'
    exit
fi

num=`wc -l $1|awk '{print $1}'`
one=`echo $num/10|bc`

sort -k7gr $1 >sorted.bed
split sorted.bed -l $one split
ls split*|xargs -n1 sh -c 'sort -k8gr $0 >$0_sorted'
cat split*sorted
rm split*

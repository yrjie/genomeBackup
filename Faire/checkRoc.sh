if [ $# -lt 2 ]
then
    echo 'Usage: bedfile TFfile'
    exit
fi

one=100000
for ((i=1;i<=10;i++))
do
    now=`echo $one*$i|bc`
    echo $now `sort -k7gr $1 |head -n $now|windowBed -a $2 -b stdin -u -w 0|wc -l`
done

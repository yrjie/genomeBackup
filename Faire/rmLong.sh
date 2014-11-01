if [ $# -lt 2 ]
then
    echo 'Usage: bedfile len'
    exit
fi

awk -v ma=$2 'BEGIN{OFS="\t"}{if ($4<ma&&$4>-ma) print $1,$2,$3,".",".",$5}' $1

if [ $# -lt 1 ]
then
    echo 'Usage: bedpe'
    exit
fi

bed=${1%%.*}.bed
awk '{printf $1"\t"$2"\t"$3"\t.\t.\t+\n"$4"\t"$5"\t"$6"\t.\t.\t-\n"}' $1 >$bed

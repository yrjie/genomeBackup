if [ $# -lt 2 ]
then
    echo 'Usage: sorted.bed extLen'
    exit
fi

newLen=`echo 100+2*$2|bc`
echo $newLen
pre=${1%%.*}
awk -v ext=$2 -v thre=$newLen 'BEGIN{OFS="\t"} $6=="+" {$3+=ext} $6=="-"&&$2>ext {$2-=ext} $5<0&&$5>-thre {next} $5>0&&$5<thre {$3=$2+$5} {print $0}' $1 >$pre"_ext"$2.bed

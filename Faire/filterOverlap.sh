if [ $# -lt 2 ]
then
    echo 'Usage: sorted.bed fraglen'
    exit
fi

pre=${1%%.*}
awk -v frag=$2 'BEGIN{OFS="\t"} $5>-frag&&$5<0 {next;} $5>0&&$5<frag {$3=$2+$5;} {print $0}' $1 >$pre"_filt"$2.bed

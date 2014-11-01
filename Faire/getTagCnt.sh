if [ $# -lt 2 ]
then
    echo 'Usage: active repress'
    exit
fi

smp=200000
rm temp.dat
cut -f1-4 $1 |uniq -c|shuf -n$smp |awk 'BEGIN{printf "Active"}{printf "\t"$1} END{printf "\n"}' >>temp.dat
cut -f1-4 $2 |uniq -c|shuf -n$smp |awk 'BEGIN{printf "Repressive"}{printf "\t"$1} END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R


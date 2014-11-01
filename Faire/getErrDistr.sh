if [ $# -lt 1 ]
then
    echo 'Usage: resultFile'
    exit
fi

awk 'BEGIN{printf "err"}{if (NR%2==0) printf "\t"$1}' $1 >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

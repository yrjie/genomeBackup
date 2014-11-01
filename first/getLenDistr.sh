if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

tmpFile=tmp$RANDOM
awk '{if (NR==1) printf "len"; printf "\t"$3-$2}' $1 >$tmpFile
R --no-save --slave --args $tmpFile <~/bin/plotBox.R
rm $tmpFile

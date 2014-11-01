if [ $# -lt 2 ]
then
    echo 'Usage: peakfile frag.bed'
    echo 'The 5th column of frag should be the length'
    exit
fi

colN=`awk 'NR==1 {print NF+5}' $1`
windowBed -a $1 -b $2 -w 0|cut -f$colN |awk '$1<0 {$1=-$1} {print $1}' >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotHist.R

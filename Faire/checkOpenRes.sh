if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    echo 'The 4th col is score, 5th col is RNA value'
    exit
fi

cut -f4,5 $1 >temp.dat
#R --no-save --slave --args temp.dat <~/bin/plotXY.R
python ~/bin/getRankCorr.py temp.dat


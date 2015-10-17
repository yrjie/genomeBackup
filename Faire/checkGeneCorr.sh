if [ $# -lt 2 ]
then
    echo 'Usage: peakfile tssSig'
    echo 'The 4th col of peakfile is the score, the last col of tssSig is the exp'
    exit
fi

tmpfile=tmp$RANDOM
#windowBed -a $1 -b $2 -w 0 |awk 'BEGIN{OFS="\t"}{print $4,$NF}' >$tmpfile
#windowBed -a $1 -b $2 -w 0 |awk 'BEGIN{OFS="\t"}{print $4,log($NF+0.01)}' >$tmpfile
#windowBed -a $1 -b $2 -w 0 |awk 'BEGIN{OFS="\t"}{print log($4+0.01),$NF}' >$tmpfile
#windowBed -a $1 -b $2 -w 0 |awk 'BEGIN{OFS="\t"}{print log($4+0.01),log($NF+0.01)}' >$tmpfile
#R --no-save --slave --args $tmpfile <~/bin/plotXY.R

windowBed -a $1 -b $2 -w 0 |awk 'BEGIN{OFS="\t"}{print $4,$NF}' >$tmpfile
python ~/bin/getRankCorr.py $tmpfile

rm $tmpfile


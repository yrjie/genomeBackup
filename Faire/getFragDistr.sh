if [ $# -lt 1 ]
then
    echo 'Usage: AtacTagFile'
    exit
fi

# ../ATAC/alignedRep3_noM_sorted.bed
# the 5th column of aligned is fraglen

tmpFile=tmp$RANDOM
colF=`awk 'NR==1{print NF-1}' $1`

#awk -v colF=$colF '{print -$colF}' $1 >$tmpFile
#awk -v colF=$colF -v lim=1000 '$colF<lim&&$colF>-lim {print -$colF}' $1 >$tmpFile
awk -v colF=$colF -v lim=1000 '$colF<0 {$colF=-$colF} {print $colF}' $1 >$tmpFile
R --no-save --slave --args $tmpFile <~/bin/plotHist.R

rm $tmpFile

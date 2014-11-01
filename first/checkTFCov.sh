if [ $# -lt 2 ]
then
    echo "Usage: sh checkTFCov.sh prefix TFfile"
    exit
fi

# $1 prefix
# $2 TF file
win=0
#wc -l $1_c*.bed
#windowBed -a $1"_c1.bed" -b $2 -u -w $win|wc -l
#windowBed -a $1"_c2.bed" -b $2 -u -w $win|wc -l
#windowBed -a $1"_c3.bed" -b $2 -u -w $win|wc -l
#windowBed -a $1"_c4.bed" -b $2 -u -w $win|wc -l
#windowBed -a $1"_c5.bed" -b $2 -u -w $win|wc -l

wc -l $1*.clust.bed
windowBed -a $1"0.clust.bed" -b $2 -u -w $win|wc -l
windowBed -a $1"1.clust.bed" -b $2 -u -w $win|wc -l
windowBed -a $1"2.clust.bed" -b $2 -u -w $win|wc -l
windowBed -a $1"3.clust.bed" -b $2 -u -w $win|wc -l
windowBed -a $1"4.clust.bed" -b $2 -u -w $win|wc -l

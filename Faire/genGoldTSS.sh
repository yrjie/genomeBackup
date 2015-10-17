if [ $# -lt 3 ]
then
    echo 'Usage: RNAFile cageFile ratio'
    echo '  The 4th column should be the exp value'
    exit
fi

rowN=`wc -l $1|awk -v ratio=$3 '{print int($1*ratio)}'`
echo $rowN
sort -k4gr $1 |head -n$rowN >temp1.bed
sort -k4gr $2 |head -n$rowN >temp2.bed

# 4 RNA, 5 CAGE, 6 ATAC
windowBed -a temp1.bed -b temp2.bed -w 0 |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$5}' >temp.bed


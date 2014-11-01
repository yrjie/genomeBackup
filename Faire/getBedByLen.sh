if [ $# -lt 1 ]
then
    echo 'Usage: aligned.bed'
    exit
fi
awk 'BEGIN{OFS="\t"}{if ($4>-50&&$4<50) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_50.bed
awk 'BEGIN{OFS="\t"}{if ($4>-100&&$4<100) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_100.bed
awk 'BEGIN{OFS="\t"}{if ($4>-150&&$4<150) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_150.bed
awk 'BEGIN{OFS="\t"}{if ($4>-200&&$4<200) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_200.bed
#awk 'BEGIN{OFS="\t"}{if ($4>-50&&$4<50) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_50.bed
#awk 'BEGIN{OFS="\t"}{if ($4>-50&&$4<50) print $1,$2,$3,".",$4,$5}' $1 >alignedRep3_50.bed

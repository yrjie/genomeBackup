if [ $3 -lt 1 ]
then
    echo 'Usage: vimIt'
    exit
fi

awk 'BEGIN{OFS="\t"}{if ($5<0) $5=-$5; if ($5>1000||$5<200) next; print $0}' alignedRep3.bed |sortBed >alignedRep3_200+HT.bed
awk 'BEGIN{OFS="\t"}{if ($5>1000||$5<200) next; print $0}' alignedRep3_fullLen.bed >alignedRep3_200+.bed

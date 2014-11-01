if [ $# -lt 3 ]
then
    echo 'Usage: pos.bed neg.bed bwfile'
    exit
fi

tmpp1=tempp1$RANDOM
tmpp2=tempp2$RANDOM
tmpn1=tempn1$RANDOM
tmpn2=tempn2$RANDOM
bigWigSummaryBatch $3 $1 1 >$tmpp1
awk '{print $3-$2}' $1 >$tmpp2
bigWigSummaryBatch $3 $2 1 >$tmpn1
awk '{print $3-$2}' $2 >$tmpn2
paste $tmpp1 $tmpp2 arff/strand_pos.tab >temp1
paste $tmpn1 $tmpn2 arff/strand_neg.tab >temp2
tab2arff.py temp1 temp2 FaireRepro
#tab2arff.py $tmpp1 $tmpn2 FaireRepro
rm $tmpp1 $tmpp2 $tmpn1 $tmpn2

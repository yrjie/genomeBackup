if [ $# -lt 2 ]
then
    echo 'Usage: peakFile colForSig'
    exit
fi

shuf -n200000 $1 >temp.bed
bigWigSummaryBatch ~/Faire/ATAC/AtacGm12878_pileup.bigWig temp.bed 1 >temp1.dat
cut -f$2 temp.bed >temp2.dat
R --no-save --slave --args temp1.dat temp2.dat <~/bin/plotXY.R

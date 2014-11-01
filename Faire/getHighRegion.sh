if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

echo $1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 >temp1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 >temp2
paste $1 temp1 temp2 |awk 'BEGIN{OFS="\t"}{a=$7+0.01;b=$8+0.01;$9=log(a/b);print $1,$2,$3,$4,$5,$6,$7,$8,$9}' |sort -k9gr >${1%%.*}_sig.bed

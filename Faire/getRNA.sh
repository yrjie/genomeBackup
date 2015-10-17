if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

echo $1
colN=`awk 'NR==1 {print NF}' $1`
file=`basename $1`
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 >temp1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 >temp2
paste $1 temp1 temp2 |awk -v colN=$colN 'BEGIN{OFS="\t"}{a=$(colN+1)+0.01;b=$(colN+2)+0.01;$(colN+3)=log(a/b);print $0}' >${file%%.*}_sig.bed

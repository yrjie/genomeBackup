if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    echo 'The 4th column should be the score'
    exit
fi

echo $1
colN=`awk 'NR==1 {print NF}' $1`
file=`basename $1`

#/gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapPlusRawRep1.bigWig
#/gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapMinusRawRep1.bigWig

/home/ruijie/gbdb/hg19/bbi/wgEncodeRikenCageGm12878NucleolusTotalMinusSignal.bigWig
/home/ruijie/gbdb/hg19/bbi/wgEncodeRikenCageGm12878NucleolusTotalPlusSignal.bigWig


bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 >temp1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 >temp2
#paste $1 temp1 temp2 |awk -v colN=$colN 'BEGIN{OFS="\t"}{a=$(colN+1)+0.01;b=$(colN+2)+0.01;tmp=log(a/b);if (tmp<0) tmp=-tmp; print $4"\t"tmp}' >temp.dat
#paste $1 temp1 temp2 |awk -v colN=$colN 'BEGIN{OFS="\t"}{a=$(colN+1)+0.01;b=$(colN+2)+0.01;tmp=log(a/b);if (tmp<0) tmp=-tmp; if (tmp>0) print $4"\t"tmp}' >temp.dat

paste $1 temp1 temp2 |awk 'BEGIN{OFS="\t"}{a=$NF+$(NF-1); print $4,a}' >temp.dat

R --no-save --slave --args temp.dat <~/bin/plotXY.R

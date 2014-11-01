if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig $1 1 >temp1
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig $1 1 >temp2
paste $1 temp1 temp2|awk 'BEGIN{OFS="\t"}{print $3-$2,$11+$12}'|sort -k1g|awk 'BEGIN{now=0}{if ($1>=now){if (now>0) printf "\n"; printf now; if (now<1000) now+=200;else if (now<10000) now+=1000;} printf "\t"log($2+0.1)}' >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

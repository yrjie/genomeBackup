if [ $# -lt 2 ]
then
    echo "bedfile cellline"
    exit
fi
echo $1
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeq$2R2x75Th1014Il200SigRep1V4.bigWig $1 5 |awk 'BEGIN {l=0;r=0}{a=$1+$2;b=$4+$5;if (a>2*b) l++; else if (b>2*a) r++}END {print l,r}'
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeq$2R2x75Th1014Il200SigRep1V4.bigWig $1 2 |awk 'BEGIN {l=0;r=0}{if ($1>2*$2) l++;else if($2>2*$1)r++;}END {print l,r}'

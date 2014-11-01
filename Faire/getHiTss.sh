if [ $# -lt 1 ]
then
    echo "Usage: tss.bed"
    exit
fi

bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig $1 1 >temp.sig
paste $1 temp.sig|sort -k7gr|head -n10000|sortBed|mergeBed >temp.bed

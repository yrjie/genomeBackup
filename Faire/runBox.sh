if [ $# -lt 4 ]
then
    echo 'Usage: Lbed Mbed Rbed cell'
    exit
fi

#mfile=/gbdb/hg19/bbi/wgEncodeGisRnaSeq$4CytosolPapMinusRawRep1.bigWig
#pfile=/gbdb/hg19/bbi/wgEncodeGisRnaSeq$4CytosolPapPlusRawRep1.bigWig

mfile=/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeq$4CellLongnonpolyaMinusRawSigRep1.bigWig
pfile=/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeq$4CellLongnonpolyaPlusRawSigRep1.bigWig

rm temp.dat
bigWigSummaryBatch $mfile $1 1 |awk '{if (NR==1) printf "LMinus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch $pfile $1 1 |awk '{if (NR==1) printf "LPlus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch $mfile $2 1 |awk '{if (NR==1) printf "MMinus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch $pfile $2 1 |awk '{if (NR==1) printf "MPlus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch $mfile $3 1 |awk '{if (NR==1) printf "RMinus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch $pfile $3 1 |awk '{if (NR==1) printf "RPlus"; printf "\t"$1}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

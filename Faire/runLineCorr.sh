if [ $# -lt 2 ]
then
    echo 'Usage: geneExp pathWay'
    exit
fi

#bw1=/home/ruijie/gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapMinusRawRep1.bigWig
#bw2=/home/ruijie/gbdb/hg19/bbi/wgEncodeGisRnaSeqGm12878CytosolPapPlusRawRep1.bigWig

bw1=/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaMinusRawSigRep1.bigWig
bw2=/gbdb/hg19/bbi/wgEncodeCshlLongRnaSeqGm12878CellLongnonpolyaPlusRawSigRep1.bigWig

#awk 'BEGIN{OFS="\t"}{a=$7+$8; print $5,a}' ../hg19Tss1000_sig.bed >geneExp.dat
#awk 'BEGIN{OFS="\t"}{a=$7+$8; if (a>100) print $5,a}' ../hg19Tss1000_sig.bed >geneExp.dat

#python /home/ruijie/Faire/open-gene/linear/getLinearCorr.py $1 /home/ruijie/Faire/open-gene/linear/ReactomePathways.gmt

python /home/ruijie/Faire/open-gene/linear/getLinearCorr.py $1 $2



#[ruijie@genome linear]$ sh runLineCorr.sh geneExp.dat posGene.xls
#4911 -0.00450880147599
#[ruijie@genome linear]$ sh runLineCorr.sh geneExp.dat negGene.xls
#774 0.227610782167

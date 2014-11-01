if [ $# -lt 2 ]
then
    echo 'Usage: Aln.bed head'
    exit
fi

head=$2
bedtools genomecov -i $1 -g ~/genome/human.hg19.genome -bg -split -strand + >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/human.hg19.genome $head.plus.bigWig
bedtools genomecov -i $1 -g ~/genome/human.hg19.genome -bg -split -strand - >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/human.hg19.genome $head.minus.bigWig
#bigWigSummaryBatch 2bp/$head.minus.bigWig data/wgEncodeAwgTfbsUtaK562CtcfUniPk_plus.narrowPeak 1 -type=max >pattern/K562Ctcf_plus.xls

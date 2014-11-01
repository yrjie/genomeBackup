if [ $# -lt 1 ]
then
    echo 'Usage: bamfile'
    exit
fi

# $1 input bam
bedF=${1%\.*}.bed
python ~/bin/extendTag.py $1 ~/genome/human.hg19.genome 150 >$bedF
num=`wc -l $bedF |cut -d' ' -f1`
echo $bedF $num
sc=`echo 1000000/$num |bc -l`
file=`basename $1`
bedtools genomecov -ibam $1 -g ~/genome/human.hg19.genome -bg -scale $sc >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/human.hg19.genome coverage/${file%\.*}.bigWig
#bigWigSummaryBatch 2bp/$head.minus.bigWig data/wgEncodeAwgTfbsUtaK562CtcfUniPk_plus.narrowPeak 1 -type=max >pattern/K562Ctcf_plus.xls

rm $bedF

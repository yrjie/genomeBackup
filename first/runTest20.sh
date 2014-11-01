if [ $# -lt 1 ]
then
    echo "Usage: FP20.bed"
    exit
fi

num=`wc -l $1 |cut -d' ' -f1`
one=`echo $num/10 |bc`
max1=10000

if [ $max1 -lt $one ]
then
    one=$max1
fi

echo "CTCF"
windowBed -a $1 -b ~/DNase/data/wgEncodeAwgTfbsUwGm12878CtcfUniPk.narrowPeak -u -w 0 >temp1.bed
sh ~/bin/getLenDistr.sh temp1.bed

echo "Cmyc"
windowBed -a $1 -b /home/ruijie/DNase/data/wgEncodeAwgTfbsUtaGm12878CmycUniPk.narrowPeak -u -w 0 >temp1.bed
sh ~/bin/getLenDistr.sh temp1.bed


#shuffleBed -i /home/browser/BASIC/import-data/broadPeak/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak -g ~/genome/human.hg19.genome >shuffled.bed
#
#sortFile=sorted$RANDOM
#sort -k7 -gr $1 >$sortFile
#
#for ((i=1;i<=10;i++)) do
#    now=`echo $one*$i|bc`
#    # be careful about the column of score
#    head -n $now $sortFile >tempFP20
#    a=`windowBed -a /home/browser/BASIC/import-data/broadPeak/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak -b tempFP20 -u -w 0|wc -l`
#    b=`windowBed -a shuffled.bed -b tempFP20 -u -w 0|wc -l`
#    echo $b $a
#done
#rm $sortFile

#windowBed -a /home/browser/BASIC/import-data/broadPeak/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak -b $1 -u -w 0|wc -l
#windowBed -a /home/browser/BASIC/import-data/broadPeak/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak -b $2 -u -w 0|wc -l

#windowBed -a shuffled.bed -b $1 -u -w 0|wc -l
#windowBed -a shuffled.bed -b $2 -u -w 0|wc -l


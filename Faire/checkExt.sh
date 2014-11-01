if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    exit
fi

echo $1

a=`windowBed -a ../ATAC/uncovTF_pile.bed -b ~/data/hg19/bed/wgEncodeBroadHmmGm12878HMM.bed -w 0|grep Strong_Enh|wc -l`

b=`windowBed -a ../ATAC/uncovTF_pile.bed -b ~/data/hg19/bed/wgEncodeBroadHmmGm12878HMM.bed -w 0|grep Strong_Enh|windowBed -a stdin -b $1 -u -w 0|wc -l`

c=`windowBed -a ../ATAC/uncovTF_pile.bed -b ~/data/hg19/bed/wgEncodeBroadHmmGm12878HMM.bed -w 0|grep Heterochro|wc -l`

d=`windowBed -a ../ATAC/uncovTF_pile.bed -b ~/data/hg19/bed/wgEncodeBroadHmmGm12878HMM.bed -w 0|grep Heterochro|windowBed -a stdin -b $1 -u -w 0|wc -l`

echo -e $a"\t"$b"\n"$c"\t"$d

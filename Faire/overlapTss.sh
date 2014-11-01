if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    exit
fi

export file=$1
wc -l $1|awk '{print $1}'
windowBed -a $1 -b ~/data/hg19/hg19Tss.bed -u -w 0|wc -l
ls high*s_tss.bed |xargs -n1 sh -c 'windowBed -a $0 -b $file -u -w 0|wc -l'
windowBed -b $1 -a ~/data/hg19/hg19Tss.bed -u -w 0|wc -l

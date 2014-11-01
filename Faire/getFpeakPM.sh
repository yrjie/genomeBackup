if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

windowBed -a $1 -b ~/data/hg19/narrowPeak/wgEncodeUwDnaseH1hescPkRep1.narrowPeak -u -w 200 |awk 'BEGIN{OFS="\t"}{for (i=11;i<=30;i++) print i-20, $i}' >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotXY.R
mv temp.png tempD.png

windowBed -a $1 -b ~/data/hg19/narrowPeak/wgEncodeUwDnaseH1hescPkRep1.narrowPeak -v -w 2000 |awk 'BEGIN{OFS="\t"}{for (i=11;i<=30;i++) print i-20, $i}' >temp.dat
R --no-save --slave --args temp.dat <~/bin/plotXY.R
mv temp.png tempnoD.png

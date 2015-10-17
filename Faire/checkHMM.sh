if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    exit
fi

cnt1=tmp1$RANDOM
cnt2=tmp2$RANDOM

shuf=shuffled$RANDOM

#hmm=../chromHMM/wgEncodeBroadHmmGm12878HMM.bed
hmm=wgEncodeBroadHmmGm12878HMM_mergeEnh.bed
colN=`awk 'NR==1 {print NF}' $1`
colSt=`echo $colN+4|bc`

shuffleBed -i $1 -g ~/genome/human.hg19.genome >$shuf

windowBed -a $1 -b $hmm -w 0 |cut -f$colSt |sort |uniq -c|awk '{print $2"\t"$1}' >$cnt1
windowBed -a $shuf -b $hmm -w 0 |cut -f$colSt |sort |uniq -c|awk '{print $2"\t"$1}' >$cnt2

paste $cnt1 $cnt2 |awk '{print $1"\t"($2+1)/($4+1)}' > temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBar.R

cut -f2 temp.dat

rm $cnt1 $cnt2 $shuf

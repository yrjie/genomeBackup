if [ $# -lt 2 ]
then
    echo 'Usage: peakFile colSig'
    exit
fi

num=`wc -l $1|awk '{print $1}'`
one=`echo $num/10|bc`

sort -k$2gr $1 |head -n$one >temp.bed
wc -l temp.bed

#ls state_*|xargs -n1 sh -c 'echo $0;windowBed -a temp.bed -b $0 -u -w 0|wc -l'
ls state_*|xargs -n1 sh -c 'echo $0;windowBed -b temp.bed -a $0 -u -w 0|wc -l'

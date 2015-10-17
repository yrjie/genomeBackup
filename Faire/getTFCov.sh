if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    exit
fi

export peak=$1
cut -f1 ../UCSCnrpk.lst |grep Gm12878|xargs -n1 sh -c 'a=`wc -l $0;` b=`windowBed -a $0 -b $peak -u -w 0|wc -l`; echo -e $a $b'|awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$3/$1}' |sort -k4gr >temp.txt

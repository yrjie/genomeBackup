if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

#cut -f1 ../UCSCnrpk.lst |grep Gm12878|xargs -n1 sh -c 'a=`wc -l $0;` b=`windowBed -a $0 -b Gm12878Enh.bed -u -w 0|wc -l`; echo -e $a $b'|awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$3/$1}' >temp.txt
cut -f1 ../UCSCnrpk.lst |grep Gm12878|xargs -n1 sh -c 'a=`wc -l $0;` b=`windowBed -b $0 -a Gm12878Enh.bed -u -w 0|wc -l`; echo -e $a $b'|awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$3/26829}' >temp1.txt

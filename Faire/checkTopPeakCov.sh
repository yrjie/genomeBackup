if [ $# -lt 2 ]
then
    echo 'Usage: peak tags_sortedLen'
    exit
fi

sortPk=sortPk$RANDOM
sort -k4gr $1 >$sortPk

pnum=`wc -l $1|cut -f1 -d' '`
pone=`echo $pnum/10|bc`
split -l$pone $sortPk psplited

#tnum=`wc -l $2|cut -f1 -d' '`
#tone=`echo $tnum/10|bc`
#split -l$tone $2 tsplited

suffix=(a b c d e f g h i j)
suffixT=(a b c d e f g h)
len=(0 100 180 247 315 473 558 615 1000000)

for ((i=0;i<8;i++))
do
    echo ${len[$i]} ${len[($i+1)]}
    awk -v a=${len[$i]} -v b=${len[($i+1)]} '$5>a&&$5<b {print $0}' $2 >tspliteda${suffixT[i]}
done

#for i in ${suffix[*]}
#do
#    echo pspliteda$i
#    for j in ${suffixT[*]}
#    do
#	windowBed -a tspliteda$j -b pspliteda$i -u -w 0|wc -l
#    done
#done

rm psplited* #tsplited*
rm $sortPk

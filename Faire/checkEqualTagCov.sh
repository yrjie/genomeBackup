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

tnum=`wc -l $2|cut -f1 -d' '`
tone=`echo $tnum/10|bc`
split -l$tone $2 tsplited

suffix=(a b c d e f g h i j)

for i in ${suffix[*]}
do
    echo pspliteda$i
    for j in ${suffix[*]}
    do
	windowBed -a tspliteda$j -b pspliteda$i -u -w 0|wc -l
    done
done

rm psplited* tsplited*
rm $sortPk

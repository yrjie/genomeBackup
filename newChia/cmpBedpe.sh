if [ $# -lt 3 ]
then
    echo 'Usage: bedpe1 bedpe2 window'
    echo 'The 7th row should be the ID'
    exit
fi

wc -l $1
wc -l $2
pairToPair -a $1 -b $2 -slop $3 -is |cut -f7|uniq|wc -l
pairToPair -b $1 -a $2 -slop $3 -is |cut -f7|uniq|wc -l
#pairToPair -b /home/ruijie/newChia/CHK159Tag_AABB/wgEncodeGisChiaPetK562Pol2InteractionsRep1.bedpe -a $1 -slop $2 -is |wc -l
#wc -l /home/ruijie/newChia/CHK159Tag_AABB/wgEncodeGisChiaPetK562Pol2InteractionsRep2.bedpe
#pairToPair -a /home/ruijie/newChia/CHK159Tag_AABB/wgEncodeGisChiaPetK562Pol2InteractionsRep2.bedpe -b $1 -slop $2 -is |wc -l
#pairToPair -b /home/ruijie/newChia/CHK159Tag_AABB/wgEncodeGisChiaPetK562Pol2InteractionsRep2.bedpe -a $1 -slop $2 -is |wc -l

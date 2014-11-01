if [ $# -lt 2 ]
then
    echo 'Usage: overlapGene.sh bedfile genefile'
    echo 'The bedfile should be 5 columns'
    exit
fi
# $1 bed file: should be 5 columns
# $2 gene file
closestBed -a $1 -b $2 -t first -d |awk 'BEGIN{OFS="\t";print "chr","start","end","name","score","geneChr","geneStart","geneEnd","geneID","geneSymbol","strand","dist2TSS","in/outGene"} {mid=int(($2+$3)/2);if ($11=="+") {dist=mid-$7} else {dist=$8-mid} if (dist>=0&&dist<=$8-$7) {inG="in"} else {inG="out"} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,dist,inG}'

if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

tmpfile='temp1'
awk 'BEGIN {OFS="\t"}{print $1,$2,$3,".","."}' $1 >$tmpfile
sh ~/bin/overlapGene.sh $tmpfile ~/data/hg19/hg19Tss.bed | awk '{if (NR==1) next;a=$(NF-1); if (a<0) a=-a; print 1/(1+a/1000.0)}'

if [ $# -lt 3 ]
then
    echo 'Usage: bed threshold suffix'
    exit
fi

# ../ATAC/alignedRep3_sorted.bed

bg=temp$3.bedGraph
bedtools genomecov -i $1 -g ~/genome/human.hg19.genome -bg >$bg
awk -v th=$2 'BEGIN{OFS="\t"}{if ($4>th) print $1,$2,$3,".",$4}' $bg |mergeBed -d 1 -scores mean >peaks/pilePeak$3.bed
#awk -v th=$2 'BEGIN{OFS="\t"}{if ($4>th) print $1,$2,$3,".",$4}' $bg |mergeBed -d 1 -scores mean >equalSplit/pilePeak$3.bed
rm $bg

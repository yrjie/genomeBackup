if [ $# -lt 2 ]
then
    echo 'Aln.bed cutoff'
    exit
fi

bedtools genomecov -i $1 -g ~/genome/human.hg19.genome -bg |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",$4}' >AtacLongFrag.bedGraph
awk -v cut=$2 '{if ($5>cut) print $0}' AtacLongFrag.bedGraph |mergeBed -d 1 -scores mean >pile.bed


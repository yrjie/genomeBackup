paste seqDepth/hg19_500.bed $1 |awk 'BEGIN {OFS="\t"} {print $1,$2,$3,".",$4}'

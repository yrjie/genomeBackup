if [ $# -lt 1 ]
then
    echo 'Usage: prefix'
    exit
fi
# $1 prefix
closestBed -a $1_c1.bed -b ~/Faire/data/knownGene.bed -t first |cut -f10|sort |uniq -c
closestBed -a $1_c2.bed -b ~/Faire/data/knownGene.bed -t first |cut -f10|sort |uniq -c
closestBed -a $1_c3.bed -b ~/Faire/data/knownGene.bed -t first |cut -f10|sort |uniq -c
closestBed -a $1_c4.bed -b ~/Faire/data/knownGene.bed -t first |cut -f10|sort |uniq -c
#closestBed -a $1_c0.bed -b data/knownGene.bed -t first |cut -f12|sort -r|uniq -c
#closestBed -a $1_c1.bed -b data/knownGene.bed -t first |cut -f12|sort -r|uniq -c
#closestBed -a $1_c2.bed -b data/knownGene.bed -t first |cut -f12|sort -r|uniq -c
#closestBed -a $1_c3.bed -b data/knownGene.bed -t first |cut -f12|sort -r|uniq -c
#closestBed -a $1_c4.bed -b data/knownGene.bed -t first |cut -f12|sort -r|uniq -c
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCaltechRnaSeqGm12878R2x75Th1014Il200SigRep1V4.bigWig data/Gm12878FaireSp1_c1.bed 2 |awk 'BEGIN {l=0;r=0}{if ($1>2*$2) l++;elif ($2>2*$1) r++}END {print l,r}'

#bigWigSummaryBatch ~/DNase/bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1.bw FaireGm12878_1no2.bed 1 >temp1
#slopBed -i FaireGm12878_1no2.bed -b 2000 -g ~/genome/human.hg19.genome >temp.bed
#
#bigWigSummaryBatch ~/DNase/bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1.bw temp.bed 1 >temp2
#
#paste temp1 temp2 |awk 'BEGIN {printf "1no2"}{printf "\t"$1/($2+0.1)}END{printf "\n"}' >temp.dat
#
#R --no-save --slave --args temp.dat "Faire Fold" <~/bin/plotBox.R

#afile='FaireGm12878_1no2.bed'
#bfile='FaireGm12878_12.bed'
afile='../ATAC/FaireGm12878_noA.bed'
bfile='../ATAC/FaireGm12878_A.bed'

bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw $afile 1 >temp1
bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw $afile 1 >temp2
paste temp1 temp2 |awk 'BEGIN {printf "Gm12878 FnoA"}{printf "\t"($1+0.1)/($2+0.1)}END{printf "\n"}' >temp.dat
bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw $bfile 1 >temp1
bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw $bfile 1 >temp2
paste temp1 temp2 |awk 'BEGIN {printf "Gm12878 FA"}{printf "\t"($1+0.1)/($2+0.1)}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat "Faire Fold" <~/bin/plotBox.R

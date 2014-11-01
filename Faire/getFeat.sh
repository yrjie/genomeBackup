mfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw"
pfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw"
#mfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.minus.bw"
#pfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.plus.bw"
#mfile="/home/ruijie/Faire/2bp/FaireK562c.minus.bw"
#pfile="/home/ruijie/Faire/2bp/FaireK562c.plus.bw"

#file='../ATAC/FaireGm12878_A.bed'
file='FaireH1hesc_D.bed'
file='FaireK562_D.bed'
file='FaireGm12878_randP.bed'
bigWigSummaryBatch ~/DNase/bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1.bw $file 1 >pos1
awk '{print $3-$2}' $file >pos2
awk 'BEGIN{OFS="\t"}{win=200;mid=($2+$3)/2;print $1,mid-win,mid+win}' $file >temp.bed
bigWigSummaryBatch $mfile temp.bed 20 >pos3
bigWigSummaryBatch $pfile temp.bed 20 >pos4

#file='../ATAC/FaireGm12878_noA.bed'
file='FaireH1hesc_noD.bed'
file='FaireK562_noD.bed'
file='FaireGm12878_randN.bed'
bigWigSummaryBatch ~/DNase/bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1.bw $file 1 >neg1
awk '{print $3-$2}' $file >neg2
awk 'BEGIN{OFS="\t"}{win=200;mid=($2+$3)/2;print $1,mid-win,mid+win}' $file >temp.bed
bigWigSummaryBatch $mfile temp.bed 20 >neg3
bigWigSummaryBatch $pfile temp.bed 20 >neg4

#paste pos1 pos2 pos3 pos4 >temp1
#paste neg1 neg2 neg3 neg4 >temp2
paste pos3 pos4 >temp1
paste neg3 neg4 >temp2

#tab2arff.py temp1 temp2 Faire_sig_len_stP >arff/len_stP_A.arff
#tab2arff.py temp1 temp2 FaireH1hesc_stP >arff/H1hesc_stP.arff
tab2arff.py temp1 temp2 FaireH1hesc_stP >arff/stP_rand.arff

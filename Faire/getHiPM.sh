awk 'BEGIN{OFS="\t"}{if ($5<$6) min=$5; else min=$6;print $1,$2,$3,$4,min}' Fpeak/test/Gm12878_noP.bed |sort -k5gr |head -n50000 >temp1.bed
sh checkFpCov.sh Gm12878 temp1.bed FP1_hiPM
paste output/Gm12878FP1_all.xls output/Gm12878FP1_hiPM_all.xls |awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$9,$9/$4,$10/$5}'|less

awk 'BEGIN{OFS="\t"}{if ($5<$6) min=$5; else min=$6;print $1,$2,$3,$4,min}' Fpeak/test/Gm12878_noP.bed |sort -k5gr |tail -n50000 >temp1.bed
sh checkFpCov.sh Gm12878 temp1.bed FP1_loPM
paste output/Gm12878FP1_all.xls output/Gm12878FP1_loPM_all.xls |awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$9,$9/$4,$10/$5}'|less

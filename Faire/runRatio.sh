sh getRatio.sh ../ATAC/FaireGm12878_A.bed FA >temp.dat
sh getRatio.sh ../ATAC/FaireGm12878_noA.bed FnoA >>temp.dat
#sh getRatio.sh /home/ruijie/Faire/refine/FaireGm12878_Fseq12.bed 12 >temp.dat
#sh getRatio.sh /home/ruijie/Faire/refine/FaireGm12878_Fseq1no2.bed 1no2 >>temp.dat
#sh getRatio.sh /home/ruijie/Faire/refine/FaireGm12878_12.bed 12 >temp.dat
#sh getRatio.sh /home/ruijie/Faire/refine/FaireGm12878_1no2.bed 1no2 >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotViolin.R

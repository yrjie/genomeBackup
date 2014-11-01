bedtools bamtobed -i bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1.bam >'temp_bam.bed'
python chopHead.py 'temp_bam.bed'
grep -P "\t\+$" >../Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bed
grep -P "\t-$" >../Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bed

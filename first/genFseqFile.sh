bedtools bamtobed -i bam/Faire/wgEncodeOpenChromFaireH1hescAlnRep1.bam >'temp_bam.bed'
python chopHead.py 'temp_bam.bed'
grep -P "\t-$" out.bed >~/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.minus.bed &
grep -P "\t\+$" out.bed > ~/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.plus.bed &

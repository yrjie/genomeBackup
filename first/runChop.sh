cell=$1
#echo 312$cell'bam'
#bedtools bamtobed -i wgEncodeUwDnase$cell'AlnRep1.bam' >$cell'bam.bed'
bedtools bamtobed -i 'bam/Nsome/'$cell'.bam' >'temp_bam.bed'
python chopHead.py 'temp_bam.bed'
bedtools bedtobam -i out.bed -g ~/genome/human.hg19.genome >'tempc.bam'
samtools sort 'tempc.bam' 'tempc_sorted'
#samtools index $cell'c_sorted.bam'
#ln -s ../DNase/$cell'c_sorted.bam' ../public_html/
#ln -s ../DNase/$cell'c_sorted.bam.bai' ../public_html/
#rm $cell'c.bam'

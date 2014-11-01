bedtools bedtobam -i GSM935897_shc_1_es_aligned.bed -g ~/genome/mouse.mm9.genome >non-targeting.bam
samtools sort non-targeting.bam non-targeting_sorted
~/bin/cufflinks-2.1.1.Linux_x86_64/cufflinks non-targeting_sorted.bam -o output
awk 'BEGIN {FS=":|\t|-";OFS="\t"} {if (NR>1){print $11,$12,$13,$18}}' output/genes.fpkm_tracking >non-targeting_exp.bed


bedtools bedtobam -i GSM935899_ell3_1_es_aligned.bed -g ~/genome/mouse.mm9.genome >ell3_knock.bam
samtools sort ell3_knock.bam ell3_knock_sorted
~/bin/cufflinks-2.1.1.Linux_x86_64/cufflinks ell3_knock_sorted.bam -o ell3_knock
awk 'BEGIN {FS=":|\t|-";OFS="\t"} {if (NR>1){print $11,$12,$13,$18}}' ell3_knock/genes.fpkm_tracking >ell3_knock_exp.bed

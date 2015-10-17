if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

Gm12878TssPro_gene.bed Gm12878TssPro_cageSum.bed |awk '$6=="+"{print $2==$8}$6=="-"{print $3==$9}' |grep 0
paste Gm12878TssPro_geneCage.bed Gm12878TssPro_cageSum.bed |awk 'BEGIN{OFS="\t"}{print $7+$8,$NF}' >temp.dat



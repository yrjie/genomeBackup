awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$8,$9,$5}' mRNA_H3k4me3StdPk.xls >mRNA_H3k4me3StdPk.bedpe
pairToPair -a mRNA_H3k4me1StdPk.bedpe -b 00min/00min.bedpe -slop 200 -is |cut -f1-7 >00min/mRNA_H3k4me1.bedpe
ls *.bedpe|xargs -n1 sh -c 'pairToPair -a $0 -b 30min/30min.bedpe  -is |cut -f1-7 >30min/$0'

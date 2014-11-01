sort -k7gr ../Fpeak/test/Gm12878_t4.bed |head -n10000 >temp.bed
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig temp.bed 1 |awk '{if (NR==1) printf "FairePeak"; printf "\t"$1}END{printf "\n"}' >temp.dat
shuffleBed -i temp.bed -g ~/genome/human.hg19.genome >tempShuf.bed
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig tempShuf.bed 1 |awk '{if (NR==1) printf "FaireRand"; printf "\t"$1}END{printf "\n"}' >>temp.dat
sort -k7gr ../Fpeak/test/AtacGm12878_all.bed |head -n10000 >temp.bed
bigWigSummaryBatch ~/public_html/data/AtacGm12878.bigWig temp.bed 1 |awk '{if (NR==1) printf "AtacPeak"; printf "\t"$1}END{printf "\n"}' >>temp.dat
shuffleBed -i temp.bed -g ~/genome/human.hg19.genome >tempShuf.bed
bigWigSummaryBatch ~/public_html/data/AtacGm12878.bigWig tempShuf.bed 1 |awk '{if (NR==1) printf "AtacRand"; printf "\t"$1}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

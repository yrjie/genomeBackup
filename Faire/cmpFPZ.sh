#bigWigSummaryBatch ~/public_html/data/AtacGm12878.bigWig ../refine/AtacGm12878ZnoFP.bed 1 |awk '{if (NR==1) printf "ZnoFP"; printf "\t"$1}END{printf "\n"}' >temp.dat
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw ~/Faire/refine/AtacGm12878ZnoFP.bed 1 |awk '{if (NR==1) printf "ZnoFP"; printf "\t"$1}END{printf "\n"}' >temp.dat
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw ~/Faire/refine/AtacGm12878FPnoZ.bed 1 |awk '{if (NR==1) printf "FPnoZ"; printf "\t"$1}END{printf "\n"}' >>temp.dat
bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw ~/Faire/refine/AtacGm12878FPnoZ_short.bed 1 |awk '{if (NR==1) printf "FPnoZ_short"; printf "\t"$1}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

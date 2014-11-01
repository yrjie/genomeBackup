binfile=hg19_bin.bed
python binChr.py ~/genome/human.hg19.genome 10000 >$binfile
bigWigSummaryBatch ../2bp/wgEncodeOpenChromFaireH1hescAlnRep1.plus.bw $binfile 1 >H1hescPlus.sig
bigWigSummaryBatch ../2bp/wgEncodeOpenChromFaireH1hescAlnRep1.minus.bw $binfile 1 >H1hescMinus.sig
paste H1hescPlus.sig H1hescMinus.sig >temp
R --no-save --slave --args temp <~/bin/plotXY.R

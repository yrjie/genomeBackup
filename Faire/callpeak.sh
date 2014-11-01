~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -p ../iff_GM12878/ $1
#~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -p ../iff_K562/ $1
cat *.wig >temp
wigToBigWig temp ~/genome/human.hg19.genome out.bw
rm *.wig
#fseq -l 800 -v -b <bff files> -p <iff files> aligments.bed

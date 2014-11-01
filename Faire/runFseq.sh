if [ $# -lt 2 ]
then
    echo 'Usage: bedfile prefix'
    exit
fi

#~/bin/fseq/bin/fseq -l 800 -v -b bff_35/ -p iff_K562/ -o Fseq/ -f 0 2bp/FaireK562c.minus.bed

#~/bin/fseq/bin/fseq -l 800 -v -b ./bff_35/ -p ./iff_K562/ ~/DNase/bam/Faire/wgEncodeOpenChromFaireK562AlnRep1_filter.bed -of npf -t 1 -o Fseq/
#~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -p ../iff_GM12878/ ~/DNase/bam/Faire/wgEncodeOpenChromFaireGm12878AlnRep1_filter.bed -of npf -t 1 -o .
#cat *.npf >FseqGm12878.bed
#rm *.npf
#~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ ~/DNase/bam/Faire/wgEncodeOpenChromFaireH1hescAlnRep1_filter.bed -of npf -t 1
#cat *.npf >FseqH1hesc.bed
#rm *.npf
#~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -p ../iff_GM12878/ ../ATAC/AtacGm12878/AtacGm12878Aln_filter.bed -of npf -t 4 -o . -f 0
~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -of npf $1
#~/bin/fseq/bin/fseq -l 800 -v -b ../bff_35/ -f 0 -of npf $1
#~/bin/fseq/bin/fseq -l 800 -v -of npf $1
cat *.npf >$2.bed
rm *.npf

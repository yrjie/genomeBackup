windowBed -a $1 -b enhancer/wgEncodeLicrHistoneEse14H3k04me1ME0129olaStdPk.broadPeak -w 0 -v >temp1
windowBed -a temp1 -b enhancer/wgEncodeLicrHistoneEse14H3k27acME0129olaStdPk.broadPeak -w 0 -v >temp2
windowBed -a temp2 -b enhancer/wgEncodeLicrTfbsEsb4P300ME0C57bl6StdPk.broadPeak -w 0 -v >temp3
windowBed -a temp3 -b UCSC_mm9_promoters_2500bp_2500bp.txt -w 0 -v >temp_filtered.bed

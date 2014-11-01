head=$1
bedtools genomecov -ibam 'tempc_sorted.bam' -g ~/genome/human.hg19.genome -bg -split -strand + >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/human.hg19.genome 2bp/Nsome/$head.plus.bw
bedtools genomecov -ibam 'tempc_sorted.bam' -g ~/genome/human.hg19.genome -bg -split -strand - >'temp.bedGraph'
bedGraphToBigWig 'temp.bedGraph' ~/genome/human.hg19.genome 2bp/Nsome/$head.minus.bw
#bigWigSummaryBatch 2bp/$head.minus.bigWig data/wgEncodeAwgTfbsUtaK562CtcfUniPk_plus.narrowPeak 1 -type=max >pattern/K562Ctcf_plus.xls

./homer/bin/makeTagDirectory peak/ ~/DNase/bam/Faire/wgEncodeOpenChromFaireH1hescAlnRep1.bam
./homer/bin/findPeaks peak/ -style groseq -o auto -tssSize 300 -minBodySize 300 -maxBodySize 1000 -tssFold 2 -bodyFold 1.5 -endFold 4

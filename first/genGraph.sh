bgfl=`basename $1`Graph
sortBed -i $1 |genomeCoverageBed -i stdin  -bg -g genome_length_hg19.txt >  bedgraph/$bgfl 

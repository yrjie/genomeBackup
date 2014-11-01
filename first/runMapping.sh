if [ $# -lt 1 ]
then
    echo 'Usage: fasta'
    exit
fi
file=`basename $1`
prefix=${file%%.*}
bwa index $1
bwa aln -l 42 -o 0  $1 ~/newChia/linker.fasta > tmp.sai
bwa samse  -n 10000000000 $1 tmp.sai linker.fasta > tmp.sam
samtools view -bS tmp.sam > tmp.bam
#samtools view tmp.bam |less

bowtie-build $1 $prefix"bwt.fasta"

bowtie -v2 -a -f -S $prefix"bwt.fasta" linker.fasta | samtools view -bS - > tmp.bam
#samtools view tmp.sam|less

#bwa index ~/newChia/CHK159_1.fasta
#bwa aln -l 42 -o 0  ~/newChia/CHK159_1.fasta ~/newChia/linker.fasta > tmp.sai
#bwa samse  -n 10000000000 CHK159_1.fasta tmp.sai linker.fasta > tmp.sam
#samtools view -bS tmp.sam > tmp.bam
#samtools view tmp.bam |less
#
#bowtie-build CHK159_1.fasta CHK159_1bwt.fasta
#
#bowtie -v2 -a -f -S CHK159_1bwt.fasta linker.fasta | samtools view -bS - > tmp.bam
#samtools view tmp.sam|less

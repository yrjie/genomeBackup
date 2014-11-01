if [ $# -lt 1 ]
then
    echo 'Usage: prefix'
    exit
fi

prefix=$1

mkdir $prefix

echo '1. unzip the fastq files'
rm /tmp/$prefix"_R1".fastq /tmp/$prefix"_R2".fastq

gunzip /home/guest/$prefix/$prefix*R1_*.gz -c >>/tmp/$prefix"_R1".fastq &
gunzip /home/guest/$prefix/$prefix*R2_*.gz -c >>/tmp/$prefix"_R2".fastq 

paste /tmp/$prefix"_R1".fastq /tmp/$prefix"_R2".fastq >/tmp/$prefix.seq

echo '# of pair-end reads'
wc -l /tmp/$prefix.seq

echo '2. combine pair-end reads to form fragments'
python ~/chiapet-pipeline-r261/script/combineFQ4.py /tmp/$prefix.seq 15 >/tmp/$prefix"combined.fastq"
awk 'NR%4==2 {print length($0)}' /tmp/$prefix"combined.fastq" >$prefix/fragLen.dat

echo '3. chop out the tags after filtering out the linkers (only AA, BB are processed)'
python ~/chiapet-pipeline-r261/script/getTag12.py /tmp/$prefix"combined.fastq" >/tmp/$prefix"Tag.fastq"
awk 'NR%4==2 {print length($0)}' /tmp/$prefix"Tag.fastq" >$prefix/tagLen.dat

echo '4. map the tags to the genome hg19 using bowtie'
ls /tmp/$prefix"Tag.fastq" >$prefix.cfg
sh ~/chiapet-pipeline-r261/script/callPeakhg19.sh $prefix.cfg

samtools view $prefix/*COMBINED.bam |LC_ALL=C sort -k1 >$prefix/$prefix"Tag.sam"
python ~/chiapet-pipeline-r261/script/sam2mapBig.py $prefix/$prefix"Tag.sam" > $prefix/$prefix.map

echo '5. run the chiapet pipeline'
python ~/chiapet-pipeline-r261/src/python/main/chiapet.py --asm hg19 --target POLII --lib $prefix --database chiapetdb --group-id $prefix --run 1-8 $prefix/$prefix".map"

rm /tmp/$prefix*


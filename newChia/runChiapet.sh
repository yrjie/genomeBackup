if [ $# -lt 1 ]
then
    echo 'Usage: prefix'
    exit
fi

prefix=$1

#mkdir $prefix
#gunzip /home/guest/$prefix/$prefix"R1_"*.gz -c >>/tmp/$prefix"_R1".fastq &
#gunzip /home/guest/$prefix/$prefix"R2_"*.gz -c >>/tmp/$prefix"_R2".fastq 
#
#paste /tmp/$prefix"_R1".fastq /tmp/$prefix"_R2".fastq >/tmp/$prefix.seq
#wc -l /tmp/$prefix.seq
#python combineFQ4.py /tmp/$prefix.seq 15 >/tmp/$prefix"combined.fastq"
#awk 'NR%4==2 {print length($0)}' /tmp/$prefix"combined.fastq" >$prefix/fragLen.dat
#python getTag12.py /tmp/$prefix"combined.fastq" >/tmp/$prefix"Tag.fastq"
#awk 'NR%4==2 {print length($0)}' /tmp/$prefix"Tag.fastq" >$prefix/tagLen.dat
#
#ls /tmp/$prefix"Tag.fastq" >$prefix.cfg
sh callPeakhg19.sh $prefix.cfg

samtools view $prefix/*COMBINED.bam |LC_ALL=C sort -k1 >$prefix/$prefix"Tag.sam"
python sam2mapBig.py $prefix/$prefix"Tag.sam" > $prefix/$prefix.map

python ~/chiapet-pipeline-r261/src/python/main/chiapet.py --asm hg19 --target POLII --lib $prefix --database chiapetdb --group-id $prefix --run 1-8 $prefix/$prefix".map"

#rm /tmp/$prefix*


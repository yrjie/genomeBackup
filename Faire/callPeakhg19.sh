if [ $# -lt 1 ]
then
    echo "Usage: cfgfile"
    exit
fi

echo ${1%%\.*}/
#python /home/sokemay/basespace/basespace-ui/basespace-ui/peakAnalyzer/jobserver/JQpeakCalling.py $1 /opt/bowtie2-2.1.0/bowtie2 ~/index/mm9/mm9 /home/sokemay/ChIPseqPipeline/LGL/genome_length_mm9.txt ${1%%\.*}/
python /home/sokemay/basespace/basespace-ui/basespace-ui/peakAnalyzer/jobserver/JQpeakCalling.py $1 /opt/bowtie2-2.1.0/bowtie2 ~/index/hg19-bowtie2/hg19 /home/sokemay/ChIPseqPipeline/LGL/genome_length_hg19.txt ${1%%\.*}/ mapOnly

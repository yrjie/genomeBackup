if [ $# -lt 2 ]
then
    echo 'Usage: peakfile tfcount'
    echo 'The first 4 columns will be kept'
    exit
fi

windowBed -a $1 -b ../../allGm12878TFrmdup.bed -w 0|cut -f1-4|uniq -c|awk -v tf=$2 'BEGIN{OFS="\t"} $1>tf {print $2,$3,$4,$5}' 

windowBed -a $1 -b ../../allGm12878TFrmdup.bed -v -w 0 

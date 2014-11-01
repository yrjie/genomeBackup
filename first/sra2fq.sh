if [ $# -lt 1 ]
then
    echo 'Usage: sraFile'
    exit
fi
file=`basename $1`
name=${file%%.*}

fastq-dump -A $name $1 -O /tmp --split-3

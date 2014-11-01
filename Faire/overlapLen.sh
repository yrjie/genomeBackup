if [ $# -lt 1 ]
then
    echo "Usage: bedfile"
    exit
fi

export file=$1
wc -l $1
ls ../peakLen/temp*FP.bed|xargs -n1 sh -c 'echo $0;windowBed -a $file -b $0 -u -w 200 |wc -l'
ls ../peakLen/temp*Z.bed|xargs -n1 sh -c 'echo $0;windowBed -a $file -b $0 -u -w 200 |wc -l'
ls ../peakLen/temp*fseq.bed|xargs -n1 sh -c 'echo $0;windowBed -a $file -b $0 -u -w 200 |wc -l'

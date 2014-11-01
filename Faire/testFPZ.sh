if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi

awk '{if ($1=="chr20") print int(($3+$2)/2)}' $1 |xargs -n1 sh runFpeak.sh -t -0.6 -f 300 ../ATAC/AtacGm12878/AtacGm12878Aln_filter.bed -pos 

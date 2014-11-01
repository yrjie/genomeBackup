if [ $# -lt 1 ]
then
    echo 'Usage: chiaClust'
    exit
fi

awk 'BEGIN{OFS="\t"}{gsub(":","\t",$0);gsub("-","\t",$0); print $1,$2,$3,$6,$7,$8,NR}' $1

if [ $# -lt 1 ]
then
    echo 'Usage: UCSCchia'
    exit
fi

awk 'BEGIN{OFS="\t"}{gsub(":","\t",$4);gsub("\\.\\.","\t",$4);gsub("-","\t",$4);gsub(",.*","",$4); print $4,NR}' $1

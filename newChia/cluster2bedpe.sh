if [ $# -lt 1 ]
then
    echo 'Usage: clusterFile'
    exit
fi

#awk 'BEGIN{OFS="\t"}{gsub(":","\t",$1);gsub("-","\t",$1);gsub(":","\t",$4);gsub("-","\t",$4);print $1,$4}' $1

awk 'BEGIN{OFS="\t"}{gsub(":","\t",$1);gsub("-","\t",$1);gsub(":","\t",$4);gsub("-","\t",$4);print $1,$4,$7,$8,$9}' $1

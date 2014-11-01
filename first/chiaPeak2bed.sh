if [ $# -lt 1 ]
then
    echo 'Usage: chiaPeak'
    exit
fi

awk '{gsub(":","\t",$0);gsub("-","\t",$0);print $0}' $1

if [ $# -lt 1 ]
then
    echo 'Usage: htmlfile'
    exit
fi

awk 'match($0,"ID_REF")&&match($0,"VALUE"){flag=1;next} flag&&length($0)==0{flag=0} flag {print}' $1

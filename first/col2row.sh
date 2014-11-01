if [ $# -lt 1 ]
then
    echo 'Usage: name'
    echo 'accepts only stdin'
    exit
fi

awk -v name=$1 'BEGIN {printf name}{printf "\t"$1}END {printf "\n"}'

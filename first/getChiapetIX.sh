if [ $# -lt 1 ]
then
    echo 'Usage: bedpe'
    exit
fi

cut -f1-6,8,18,19 $1

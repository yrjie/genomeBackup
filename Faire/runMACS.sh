if [ $# -lt 2 ]
then
    echo 'Usage: bedfile jobname'
    exit
fi

macs14 -t $1 -f auto -n ./$2 -p 0.005

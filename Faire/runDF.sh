if [ $# -lt 2 ]
then
    echo 'Usage: bedfile outfile'
    exit
fi

export PATH=$PATH:~/bin/DFilter/

run_dfilter.sh -d=$1 -o=$2 -std=2 -bs=100, -ks=50 

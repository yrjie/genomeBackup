if [ $# -lt 1 ]
then
    echo 'Usage: TFpeak'
    exit
fi

export TFfile=$1

ls equalSplit/pilePeaksplit*.bed |xargs -n1 sh -c 'windowBed -a $TFfile -b $0 -u -w 0|wc -l'

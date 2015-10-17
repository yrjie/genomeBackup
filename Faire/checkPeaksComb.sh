if [ $# -lt 1 ]
then
    echo 'Usage: prefix TFfile'
    exit
fi

export file=$1
export tf=$2
#ls allComb_*| xargs -n1 sh -c 'echo $0; windowBed -a $file -b $0 -u -w 0|wc -l; windowBed -a $file -b $0 -u -w 0|windowBed -a $tf -b stdin -u -w 0|wc -l; windowBed -a $file -b $0 -u -w 0|windowBed -b $tf -a stdin -u -w 0|wc -l'

ls $1* |xargs -n1 sh -c 'echo $0; wc -l $0;  windowBed -a $tf -b $0 -u -w 0|wc -l; windowBed -b $tf -a $0 -u -w 0|wc -l;'

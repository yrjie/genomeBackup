if [ $# -lt 1 ]
then
    echo 'Usage: bedfile'
    exit
fi
export file=$1
grep BroadHistoneGm12878 ../UCSCbrpk.lst |cut -f1|xargs -n1 sh -c 'echo $0;windowBed -a $file -b $0 -u -w|wc -l'

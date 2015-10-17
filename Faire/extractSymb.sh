if [ $# -lt 2 ]
then
    echo 'Usage: temp.bed tssFile'
    exit
fi

windowBed -a $1 -b $2 -w 0 |awk '$1==$7&&$2==$8&&$3==$9{print $11"\t"$4}'

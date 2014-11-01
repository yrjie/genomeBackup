if [ $# -lt 1 ]
then
    echo "Usage: Fairefile"
    exit
fi
file=`basename $1`
awk -v fname=$file 'BEGIN {split(fname,a,"Faire"); split(a[2],b,"Pk");name=b[1];OFS="\t"} {$4=name; print $0}' $1 >temp
mv temp $1

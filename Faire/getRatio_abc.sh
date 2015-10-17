if [ $# -lt 1 ]
then
    echo 'Usage: tagPkFile'
    exit
fi

echo $1
awk 'BEGIN {a=b=c=0;OFS="\t"} $5<=85 {a+=1} $5>85&&$5<=164 {b+=1} $5>164 {c+=1} END {n=a+b+c; print 1.0*a/n,1.0*b/n,1.0*c/n}' $1

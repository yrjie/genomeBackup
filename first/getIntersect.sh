chipseq=$1
dnase=$2
faire=$3
prefix=$4
echo $prefix
windowBed -a $chipseq -b $dnase -w 200 -u >intersect/temp
windowBed -a intersect/temp -b $faire -v -w 0 >intersect/$prefix'Dnase.bed'
windowBed -a intersect/temp -b $faire -u -w 0 >intersect/$prefix'DF.bed'

windowBed -a $chipseq -b $faire -w 200 -u >intersect/temp
windowBed -a intersect/temp -b $dnase -v -w 0 >intersect/$prefix'Faire.bed'

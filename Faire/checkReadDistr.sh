if [ $# -lt 2 ]
then
    echo 'Usage: peakfile readfile'
    exit
fi
num=`wc -l $1|cut -d" " -f1`
echo $num
plus=`windowBed -a $2 -b $1 -w 0 -u|grep '+$'|wc -l`
echo $plus
minus=`windowBed -a $2 -b $1 -w 0 -u|grep '\-$'|wc -l`
echo $minus
echo $plus/$num| bc -l
echo $minus/$num| bc -l

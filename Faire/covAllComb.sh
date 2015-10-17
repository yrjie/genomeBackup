if [ $# -lt 1 ]
then
    echo 'Usage: peakfile'
    echo 'show the enrichment of TF being covered'
    exit
fi

# input: split by sig
#echo -e '\n'
wc -l $1
for i in pileall3_*
do
#    echo $i
    a=`windowBed -a $1 -b $i -u -w 0|wc -l`
    #b=`windowBed -a $1 -b $i -u -w 0|windowBed -a stdin -b ../../allGm12878TFrmdup.bed -u -w 0|wc -l`
    b=`windowBed -a $1 -b $i -u -w 0|windowBed -b stdin -a ../../allGm12878TFrmdup.bed -u -w 0|wc -l`
#    echo $a
    echo $b/$a |bc -l
done

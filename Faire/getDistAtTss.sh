# output: (plus, minus)
if [ $# -lt 2 ]
then
    echo 'Usage: Tssfile swap.bedpe'
    exit
fi

windowBed -a $1 -b $2 -w 1000 | awk 'BEGIN{OFS="\t"}{midP=int(($5+$6)/2); midM=int(($8+$9)/2); print midP-$2,$10}' >temp.dat

R --no-save --slave --args temp.dat <~/bin/plotXY.R

sh ~/bin/binDataToBox.sh temp.dat 20

if [ $# -lt 1 ]
then
    echo "Usage: sh runPlot.sh prefix"
    exit
fi

num=`wc -l $1 |cut -d' ' -f1`
sc=1
if [ $num -gt 900 ]
then
    sc=`echo $num/900 |bc`
fi
cut -f2 -d" " $1.membership >temp
paste $1 temp |sort -k5|awk -v rate=$sc 'BEGIN {OFS="\t";print ".","2i/serum","281KD/NonT","Erk/WT"} {$NF="";if (NR%rate==0) print $0}' >matPlot.dat
~/bin/matrix2png -z -data matPlot.dat -mincolor blue -maxcolor green -missingcolor grey -size 30:1 -c >$1.png

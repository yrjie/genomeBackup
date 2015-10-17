if [ $# -lt 1 ]
then
    echo 'Usage: AtacTag'
    exit
fi

#python getFragRatioLR.py $1 >temp.dat
#python getFragRatioLR.py $1 |cut -f1-41 >temp.dat
python getFragRatioLR.py $1 |cut -f1-21 >temp.dat
./seq_main -i temp.dat -n 5
cut -f1-3 $1 |uniq >temp1.bed
paste temp1.bed temp.dat.membership |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' >tempC.bed


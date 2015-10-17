if [ $# -lt 2 ]
then
    echo 'Usage: datfile clustId'
    exit
fi

paste $1 $1.membership |awk -v id=$2 '$NF==id {printf NR;for (i=1;i<=NF-3;i++) printf "\t"$i; printf "\n"}' >temp1.dat
~/bin/matrix2png -z -data temp1.dat -mincolor black -maxcolor white -missingcolor grey -z -size 3:3 >temp.png


#python getFragRatioLR.py AtacTag_highPlus_tss1bp.bed |sort -k41gr |cut -f1-40 >temp.dat
#~/bin/matrix2png -z -data temp.dat -mincolor black -maxcolor white -z -missingcolor grey -size 3:3 >temp.png

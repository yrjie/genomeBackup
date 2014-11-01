if [ $# -lt 5 ]
then
    echo "Usage: file1 name1 file2 name2 bwfile"
    exit
fi
afile=$1
aname=$2
bfile=$3
bname=$4
bwfile=$5
bigWigSummaryBatch $bwfile $afile 1 |awk -v name=$aname 'BEGIN{printf name}{printf "\t"$1}END{printf "\n"}' >temp.dat
bigWigSummaryBatch $bwfile $bfile 1 |awk -v name=$bname 'BEGIN{printf name}{printf "\t"$1}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

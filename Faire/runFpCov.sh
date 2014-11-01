if [ $# -lt 3 ]
then
    echo "Usage: bed name cell"
    exit
fi

cell=$3

# for zinba and Atac
len=(0 200 400 600 1000 1500 2000 3000 4000 5000 1000000)

# for fseq
#len=(0 50 100 150 200 250 300 350 400 1000 1500 1000000)
letter=(a b c d e f g h i j k l m n o p q)

num=${#len[*]}
rm temp.dat
for ((i=1;i<$num;i++))
do
#    echo ${len[i-1]},$2_${len[i]}
    echo $2_${len[i-1]}
    awk -v a=${len[i-1]} -v b=${len[i]} '{l=$3-$2;if (l>a&&l<b) print $0}' $1 >temp.bed
    numL=`wc -l temp.bed|cut -d" " -f1`
    echo $numL
    if [ $numL -lt 1 ]
    then
	continue
    fi
    sh checkFpCov.sh $cell temp.bed $2_${letter[i-1]}${len[i-1]}
    outfile=output/$cell$2_${letter[i-1]}${len[i-1]}_all.xls
    awk -v name=${len[i-1]} 'BEGIN{printf name}{printf "\t"$1}END{printf "\n"}' $outfile >>temp.dat
    awk -v name=${len[i-1]}"Bg" 'BEGIN{printf name}{printf "\t"$2}END{printf "\n"}' $outfile >>temp.dat
done
paste output/$cell'Name.lst' output/$cell$2*_all.xls >output/$cell$2_combined.xls
rm output/$cell$2*_all.xls
R --no-save --slave --args temp.dat <~/bin/plotBox.R

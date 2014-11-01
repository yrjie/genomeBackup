if [ $# -lt 1 ]
then
    echo "Usage: fasta"
    exit
fi

headN=(pr1 pr2)
tailN=(pr1-rev pr2-rev)

lnkN=(AA AA-rev BB BB-rev)

headPr=(
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATC
)
tailPr=(
GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
)

headPr15=(
GACGCTCTTCCGATC
ACCGCTCTTCCGATC
)

tailPr15=(
GATCGGAAGAGCGTC
GATCGGAAGAGCGGT
)

lnk=(GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC
GCTGTTAAGGACGGCCAAGGCGGATGTACGGTCCTTAACAGC
GGTCAATGTCACCGTACATCCGCCTTGGCCGTGACATTGACC
GGTCAATGTCACGGCCAAGGCGGATGTACGGTGACATTGACC
)

file5=taglen5.dat
file3=taglen3.dat
rm $file5 $file3
for ((i=0;i<2;i++))
do
#    x=${headPr[$i]}
    x=${headPr15[$i]}
    for ((j=0;j<4;j++))
    do
	y=${lnk[$j]}
	echo ${headN[$i]}.*${lnkN[$j]}
	#echo $x.*$y
	#grep $x.*$y $1|awk '{print and(16,$2)}'|sort -gr|uniq -c
	grep $x.*$y $1|wc -l
	grep -Po '(?<='$x').*(?='$y')' $1|awk '{print length($0)}' >>$file5
	grep -P '(?<='$x').*(?='$y')' $1|awk '{print length($0)}' >>taglen5_frag.dat
    done
done

for ((i=0;i<2;i++))
do
    #x=${tailPr[$i]}
    x=${tailPr15[$i]}
    for ((j=0;j<4;j++))
    do
	y=${lnk[$j]}
	echo ${lnkN[$j]}.*${tailN[$i]}
	#echo $y.*$x
	#grep $y.*$x $1|awk '{print and(16,$2)}'|sort -gr|uniq -c
	grep $y.*$x $1|wc -l
	grep -Po '(?<='$y').*(?='$x')' $1|awk '{print length($0)}' >>$file3
	grep -P '(?<='$y').*(?='$x')' $1|awk '{print length($0)}' >>taglen3_frag.dat
    done
done

awk '{if (NR==1) printf "5end"; printf "\t"$0}END{printf "\n"}' $file5 >temp.dat
awk '{if (NR==1) printf "3end"; printf "\t"$0}END{printf "\n"}' $file3 >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

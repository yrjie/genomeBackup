if [ $# -lt 1 ]
then
    echo 'Usage: seqfile'
    exit
fi


lnk=(GCTGTTAAGGACCGTACATCCGCCTTGGCCGTCCTTAACAGC
GCTGTTAAGGACCGTACATCCGCCTTGGCCGTGACATTGACC
GGTCAATGTCACCGTACATCCGCCTTGGCCGTCCTTAACAGC
GGTCAATGTCACCGTACATCCGCCTTGGCCGTGACATTGACC
GCTGTTAAGGACGGCCAAGGCGGATGTACGGTCCTTAACAGC
GGTCAATGTCACGGCCAAGGCGGATGTACGGTCCTTAACAGC
GCTGTTAAGGACGGCCAAGGCGGATGTACGGTGACATTGACC
GGTCAATGTCACGGCCAAGGCGGATGTACGGTGACATTGACC
)

pr15=(GACGCTCTTCCGATC
ACCGCTCTTCCGATC
GATCGGAAGAGCGTC
GATCGGAAGAGCGGT
)

for x in ${lnk[*]}
#for x in ${pr15[*]}
do
    grep $x $1 |wc -l
done

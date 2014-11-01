if [ $# -lt 1 ]
then
    echo 'Usage: samfile'
    exit
fi

leftAdp=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC
rightAdp=CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATC
internal=CGTACATCCGCCTTGGCC
A=GCTGTTAAGGAC
Ar=GTCCTTAACAGC
B=GTGACATTGACC
Br=GGTCAATGTCAC
lst=(GCTGTTAAGGAC
GTCCTTAACAGC
GTGACATTGACC
GGTCAATGTCAC)

lnks=($A$internal$Ar $A$internal$B $Br$internal$Ar $Br$internal$B $leftAdp $rightAdp)

for x in ${lnks[*]}
do
    echo $x
    grep $x $1 |wc -l
done


#x=${lst[$1]}
#echo $x
#for y in ${lst[*]}
#do
#    echo $x$internal$y
#done
#echo $A$internal$A
#echo $A$internal$Ar
#echo $A$internal$B
#echo $A$internal$Br
#
#echo $Ar$internal$A
#echo $Ar$internal$Ar
#echo $Ar$internal$B
#echo $Ar$internal$Br
#
#echo $B$internal$A
#echo $B$internal$Ar
#echo $B$internal$B
#echo $B$internal$Br
#
#echo $Br$internal$A
#echo $Br$internal$Ar
#echo $Br$internal$B
#echo $Br$internal$Br

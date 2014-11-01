if [ $# -lt 2 ]
then
    echo "Usage: Afile Bfile"
    exit
fi
A=$1
B=$2
#A='../ATAC/FaireGm12878_A.bed'
#B='../ATAC/FaireGm12878_noA.bed'
anum=`wc -l $A|cut -f1 -d" "`
bnum=`wc -l $B|cut -f1 -d" "`
for ((i=1;i<=15;i++)) do
    awk -v now=$i '{if ($6>now) print $0}' combinedF_2.bed >temp
    #awk -v now=$i '{if ($5>now) print $0}' temp1 >temp
    #wc -l temp
    a1=`windowBed -a $A -b temp -u -w 0|wc -l`
    b1=`windowBed -a $B -b temp -u -w 0 |wc -l`
    covA=`echo $a1/$anum|bc -l`
    covB=`echo $b1/$bnum|bc -l`
    let "j=$i+1"
    echo -e $j"\t"$covA"\t"$covB
done

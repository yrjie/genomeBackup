if [ $# -lt 1 ]
then
    echo 'Usage: AlnBed'
    exit
fi
JAVAOPTS="-Xms100M -Xmx16000M"
#JAVAOPTS="-Xmx8000M"

outdir="test"
#outdir="K562"
java $JAVAOPS -cp Fpeak.jar:lib/commons-cli-1.1.jar fpeak.Main -v -b ../bff_35/ -of npf -o $outdir $*
#java $JAVAOPS -cp Fpeak.jar:lib/commons-cli-1.1.jar fpeak.Main -v -l 800 -of npf -o $outdir $*
cat $outdir/*.npf >$outdir/temp.bed
rm $outdir/*.npf
#sh runFpeak.sh ../ATAC/AtacGm12878Rep3/AtacGm12878AlnRep3.bed -t -1 -l 800 -f 300 -clip
#sh runFpeak.sh ../ATAC/AtacGm12878Rep1/AtacGm12878AlnRep1.bed -t 0.3 -l 800 -f 300 -clip
#sh runFpeak.sh ../ATAC/alignedRep3_200-300.bed -t 3 -l 100 -f 300

if [ $# -lt 2 ]
then
    echo "Usage: bedGraph cutoff"
    exit
fi

awk -v cut=$2 '{if ($5>cut) print $0}' $1 |mergeBed -d 1 -scores mean >nsome.bed
wc -l nsome.bed

windowBed -a ../Fpeak/test/AtacGm12878St1.bed -b nsome.bed -u -w 100 >temp.bed
echo 'hi with nsome'
wc -l temp.bed
windowBed -a temp.bed -b ../allGm12878TF.bed -u -w 0|wc -l

windowBed -a ../Fpeak/test/AtacGm12878St1.bed -b nsome.bed -v -w 200 >temp.bed
echo 'hi without nsome'
wc -l temp.bed
windowBed -a temp.bed -b ../allGm12878TF.bed -u -w 0|wc -l

windowBed -a ../ATAC/AtacGm12878_lowSig.bed -b nsome.bed -u -w 100 >temp.bed
echo 'low with nsome'
wc -l temp.bed
windowBed -a temp.bed -b ../allGm12878TF.bed -u -w 0|wc -l

windowBed -a ../ATAC/AtacGm12878_lowSig.bed -b nsome.bed -v -w 200 >temp.bed
echo 'low without nsome'
wc -l temp.bed
windowBed -a temp.bed -b ../allGm12878TF.bed -u -w 0|wc -l

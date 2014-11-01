if [ $# -lt 4 ]
then
    echo 'Usage: peak.lst A.bed B.bed suffix [outdir]'
    exit
fi
export peakLst=$1
export Afile=$2
export Bfile=$3
export suffix=$4
outdir='./'
if [ $# -gt 4 ]
then
    outdir=$5'/'
fi
wc -l $peakLst
wc -l $Afile
wc -l $Bfile
Anum=`wc -l $Afile|cut -d' ' -f1`
Bnum=`wc -l $Bfile|cut -d' ' -f1`
#sh -c 'echo $Afile $Bfile'
cut -f1 $peakLst|xargs -n1 sh -c 'windowBed -a $Afile -b $0 -u -w 200|wc -l' >tempA$suffix
cut -f1 $peakLst|xargs -n1 sh -c 'windowBed -a $Bfile -b $0 -u -w 200|wc -l' >tempB$suffix
paste $peakLst tempA$suffix tempB$suffix|awk -v Anum=$Anum -v Bnum=$Bnum 'BEGIN{FS="\t";OFS="\t"}{a=$3/Anum;b=$4/Bnum;print $2,a,b,a-b}'|sort -k4gr -t$'\t' >$outdir"peakDiff"$suffix.xls
rm temp*$suffix

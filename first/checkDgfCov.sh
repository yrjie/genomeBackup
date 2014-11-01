if [ $# -lt 2 ]
then
    echo "Usage: sh checkDgfCov.sh cell dnaseHotspot "
    exit
fi
cell=$1
dgf=$2
grep $cell tfbsUnipk.lst |cut -f2 >temp1
grep $cell tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'wc -l data/$0.narrowPeak'|xargs -n2 sh -c 'echo $0' >temp2
grep $cell tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'windowBed -a data/$0.narrowPeak -b '$dgf' -u -w 200|wc -l' >temp3
grep $cell tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'shuffleBed -i data/$0.narrowPeak -g ~/genome/human.hg19.genome >shuffled;windowBed -a shuffled -b '$dgf' -u -w 200|wc -l' >temp4
paste temp1 temp2 temp3 temp4 >temp
awk 'BEGIN {FS="\t"} {OFS="\t";if ($2>0) {print $1,$2,$3,$3/$2,$4/$2;} else {print $1,$2,$3,$2,$2;}}' temp >hotspot/output/$cell'FaireFP.xls'

#grep K562 tfbsUnipk.lst |cut -f1|xargs -n1 -I x1 windowBed -a data/x1.narrowPeak -b hotspot/Hotspotv.bed -u -w 200|cut -f5|awk 'BEGIN {sum=0} {sum+=$1} END {OFS="\t";print NR,sum/NR}' >temp3
#windowBed -a data/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak -b hotspot/Hotspotu.bed -w 200 -u|cut -f5|awk 'BEGIN {sum=0} {sum+=$1} END {OFS="\t";print NR,sum/NR}' >temp3

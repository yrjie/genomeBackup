if [ $# -lt 3 ]
then
    echo "Usage: sh checkFpCov.sh cell peak.bed suffix"
    exit
fi
cell=$1
dgf=$2
suffix=$3
grep $cell ~/DNase/tfbsUnipk.lst |cut -f2 >temp1$suffix
grep $cell ~/DNase/tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'wc -l ~/DNase/data/$0.narrowPeak'|xargs -n2 sh -c 'echo $0' >temp2$suffix
grep $cell ~/DNase/tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'windowBed -a ~/DNase/data/$0.narrowPeak -b '$dgf' -u -w 200|wc -l' >temp3$suffix
grep $cell ~/DNase/tfbsUnipk.lst |cut -f1|xargs -n1 sh -c 'shuffleBed -i ~/DNase/data/$0.narrowPeak -g ~/genome/human.hg19.genome >shuffled'$suffix';windowBed -a shuffled'$suffix' -b '$dgf' -u -w 200|wc -l' >temp4$suffix
#paste temp1$suffix temp2$suffix temp3$suffix temp4$suffix >output/$cell$suffix'_all'.xls
paste temp1$suffix temp2$suffix temp3$suffix temp4$suffix >temp$suffix
awk 'BEGIN {FS="\t"} {OFS="\t";if ($2>0) {print $3/$2,$4/$2;} else {print $2,$2;}}' temp$suffix >output/$cell$suffix'_all'.xls
rm temp*$suffix shuffled$suffix

#grep K562 tfbsUnipk.lst |cut -f1|xargs -n1 -I x1 windowBed -a data/x1.narrowPeak -b hotspot/Hotspotv.bed -u -w 200|cut -f5|awk 'BEGIN {sum=0} {sum+=$1} END {OFS="\t";print NR,sum/NR}' >temp3
#windowBed -a data/wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk.narrowPeak -b hotspot/Hotspotu.bed -w 200 -u|cut -f5|awk 'BEGIN {sum=0} {sum+=$1} END {OFS="\t";print NR,sum/NR}' >temp3

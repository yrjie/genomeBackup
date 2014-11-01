if [ $# -lt 2 ]
then
    echo "Usage: short.bed long.bed"
    exit
fi
cp $1 temp1.bed
cp $2 temp2.bed
windowBed -a /home/browser/BASIC/import-data/bed/wgEncodeBroadHmmGm12878HMM.bed -b temp1.bed -u -w 0|cut -f4|sort |uniq -c >temp1Chrom.tab
windowBed -a /home/browser/BASIC/import-data/bed/wgEncodeBroadHmmGm12878HMM.bed -b temp2.bed -u -w 0|cut -f4|sort |uniq -c >temp2Chrom.tab
python getDiffChrom.py temp1Chrom.tab allChrom.tab >temp1ChromR.tab
python getDiffChrom.py temp2Chrom.tab allChrom.tab >temp2ChromR.tab
R --no-save --slave --args temp1ChromR.tab <~/bin/plotBar.R
mv temp.png short.png
R --no-save --slave --args temp2ChromR.tab <~/bin/plotBar.R
mv temp.png long.png

shuffleBed -i $1 -g ~/genome/human.hg19.genome >temp1.bed
shuffleBed -i $2 -g ~/genome/human.hg19.genome >temp2.bed
windowBed -a /home/browser/BASIC/import-data/bed/wgEncodeBroadHmmGm12878HMM.bed -b temp1.bed -u -w 0|cut -f4|sort |uniq -c >temp1Chrom.tab
windowBed -a /home/browser/BASIC/import-data/bed/wgEncodeBroadHmmGm12878HMM.bed -b temp2.bed -u -w 0|cut -f4|sort |uniq -c >temp2Chrom.tab
python getDiffChrom.py temp1Chrom.tab allChrom.tab >temp1ChromR.tab
python getDiffChrom.py temp2Chrom.tab allChrom.tab >temp2ChromR.tab
R --no-save --slave --args temp1ChromR.tab <~/bin/plotBar.R
mv temp.png shortRand.png
R --no-save --slave --args temp2ChromR.tab <~/bin/plotBar.R
mv temp.png longRand.png

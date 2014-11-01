cut -f1 ~/chiapet-pipeline-r261/data/chrdata/gap/hg19.txt |uniq -c >temp2
addGapHeader.py ~/genome/human.hg19.genome temp2 ~/chiapet-pipeline-r261/data/chrdata/gap/hg19.txt >hg19.txt
mv hg19.txt ~/chiapet-pipeline-r261/data/chrdata/gap/

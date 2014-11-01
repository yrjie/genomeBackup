if [ $# -lt 2 ]
then
    echo "Usage: sh runHistone.sh prefix cell-line"
    exit
fi
#$1 prefix
#$2 cell line
python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1_c1.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1_c2.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1_c3.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1_c4.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1_c5.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
#python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1F.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/
#python /home/sokemay/ChIPseqPipeline/histonePlot/plotHistoneDist.py $1NF.bed BroadHistone$2 /home/sokemay/ChIPseqPipeline/broad_histone/hg19/ histone/

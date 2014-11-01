if [ $# -lt 1 ]
then
    echo "Usage: cell"
    exit
fi
cell=$1
sh getRoc.sh Fpeak/output/$cell'Peak.bed' $cell FP &
sh getRoc.sh MACS/Faire$cell'Allpeak.bed' $cell M &
sh getRoc.sh DFilter/track-Faire$cell'_DF.bed' $cell DF &

if [ $# -lt 2 ]
then
    echo "Usage: posBed allBed"
    exit
fi

pfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.plus.bw"
mfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireH1hescAlnRep1.minus.bw"
#pfile="/home/ruijie/Faire/2bp/FaireK562c.plus.bw"
#mfile="/home/ruijie/Faire/2bp/FaireK562c.minus.bw"
#pfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw"
#mfile="/home/ruijie/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw"

bigWigSummaryBatch $pfile $1 20 >tempp
bigWigSummaryBatch $mfile $1 20 >tempm
paste tempp tempm >pos.tab

bigWigSummaryBatch $pfile $2 20 >tempp
bigWigSummaryBatch $mfile $2 20 >tempm
paste tempp tempm >all.tab

#python getDist2Pos.py pos.tab all.tab 1 >temp1
python getDist2Pos.py pos.tab all.tab >temp1

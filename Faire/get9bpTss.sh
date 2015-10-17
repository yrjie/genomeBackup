if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

windowBed -a ../open-gene/highMinusRatio_tss.bed -b 9bp_highGene.bed -w 0|awk '{a=($2+$3)/2;b=($11+$12)/2;print b-a}' >temp.dat

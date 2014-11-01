if [ $# -lt 2 ]
then
    echo 'Usage: peakBed prefix'
    exit
fi

prefix=$2
hmm=wgEncodeBroadHmmGm12878HMM.bed

ncol=`awk 'NR==1{print NF}' $1`

windowBed -a $1 -b $hmm -w 0 | awk 'match($0, "Active_Promoter")||match($0, "4_Strong_Enhancer") {print $0}'|cut -f1-$ncol >$prefix"hmmPosPk.bed"
										windowBed -a $1 -b $hmm -w 0 | awk 'match($0, "Heterochrom")||match($0, "Insulator") {print $0}'|cut -f1-$ncol >$prefix"hmmNegPk.bed"

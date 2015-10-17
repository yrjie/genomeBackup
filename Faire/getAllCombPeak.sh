if [ $# -lt 2 ]
then
    echo 'Usage: peakTag prefix'
    exit
fi

tmpFile=tmp$RANDOM
prefix=$2

awk 'BEGIN {OFS="\t"} $4>0&&$5==0&&$6==0 {print $0,1}' >$prefix"_aUnique.bed" $1
awk 'BEGIN {OFS="\t"} $4>0&&$5>0&&$6==0 {print $0,5}' >$prefix"_ab.bed" $1
awk 'BEGIN {OFS="\t"} $4>0&&$5==0&&$6>0 {print $0,3}' >$prefix"_ac.bed" $1
awk 'BEGIN {OFS="\t"} $4>0&&$5>0&&$6>0 {print $0,7}' >$prefix"_abc.bed" $1
awk 'BEGIN {OFS="\t"} $4==0&&$5>0&&$6==0 {print $0,2}' >$prefix"_bUnique.bed" $1
awk 'BEGIN {OFS="\t"} $4==0&&$5==0&&$6>0 {print $0,4}' >$prefix"_cUnique.bed" $1
awk 'BEGIN {OFS="\t"} $4==0&&$5>0&&$6>0 {print $0,6}' >$prefix"_bc.bed" $1

#python /home/ruijie/Faire/openEM/getFragRatio.py $1 >$tmpFile
#
#awk '$4>0&&$5==0&&$6==0 {print $0}' >$prefix"_aUnique.bed" $tmpFile
#awk '$4>0&&$5>0&&$6==0 {print $0}' >$prefix"_ab.bed" $tmpFile
#awk '$4>0&&$5==0&&$6>0 {print $0}' >$prefix"_ac.bed" $tmpFile
#awk '$4>0&&$5>0&&$6>0 {print $0}' >$prefix"_abc.bed" $tmpFile
#awk '$4==0&&$5>0&&$6==0 {print $0}' >$prefix"_bUnique.bed" $tmpFile
#awk '$4==0&&$5==0&&$6>0 {print $0}' >$prefix"_cUnique.bed" $tmpFile
#awk '$4==0&&$5>0&&$6>0 {print $0}' >$prefix"_bc.bed" $tmpFile
#
#rm $tmpFile

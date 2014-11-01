if [ $# -lt 1 ]
then
    echo 'Usage: bedfile [name]'
    exit
fi

name=temp
if [ $# -gt 1 ]
then
    name=$2
fi

win=200
#cell=Gm12878
cell=H1hesc

#pfile=~/Faire/2bp/wgEncodeOpenChromFaire$cell"AlnRep1.plus.bw"
#mfile=~/Faire/2bp/wgEncodeOpenChromFaire$cell"AlnRep1.minus.bw"
pfile=~/Faire/2bp/AtacGm12878.plus.bw
mfile=~/Faire/2bp/AtacGm12878.minus.bw

#awk -v win=$win 'BEGIN{OFS="\t"}{print $1,$2-win,$3+win}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=int(($2+$3)/2);print $1,mid-win,mid+win}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{a=$2-win;b=$2+win;print $1,a,b}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=($2+$3)/2;print $1,$2-win,$2}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=($2+$3)/2;print $1,$3,$3+win}' $1 >temp2
#bigWigSummaryBatch ~/Faire/2bp/AtacGm12878.plus.bw temp1 20 >tempp
#bigWigSummaryBatch ~/Faire/2bp/AtacGm12878.minus.bw temp1 20 >tempm

#bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw temp1 20 >tempp
#bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw temp1 20 >tempm

#cp $1 temp1
#bigWigSummaryBatch $pfile temp1 20 >tempp
#bigWigSummaryBatch $mfile temp1 20 >tempm


#awk -v win=$win 'BEGIN{OFS="\t"}{a=$2-win/2;b=$2+win/2;print $1,a,b}' $1 >temp1
awk -v win=$win 'BEGIN{OFS="\t"}{a=$2-win;b=$2;print $1,a,b}' $1 >temp1
bigWigSummaryBatch $pfile temp1 20 >tempp1
bigWigSummaryBatch $mfile temp1 20 >tempm1

#awk -v win=$win 'BEGIN{OFS="\t"}{a=$3-win/2;b=$3+win/2;print $1,a,b}' $1 >temp1
awk -v win=$win 'BEGIN{OFS="\t"}{a=$3;b=$3+win;print $1,a,b}' $1 >temp1
bigWigSummaryBatch $pfile temp1 20 >tempp2
bigWigSummaryBatch $mfile temp1 20 >tempm2

paste tempp1 tempp2 >tempp
paste tempm1 tempm2 >tempm

#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig temp1 1 >temp1.sig
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig temp2 1 >temp2.sig
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig temp1 1 >temp1.sig
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig temp2 1 >temp2.sig
#bigWigSummaryBatch /gbdb/hg19/bbi/wgEncodeOpenChromFaireGm12878BaseOverlapSignal.bigWig $1 1 >temp.sig
#paste temp1.sig temp.sig temp2.sig >tempm
#bigWigSummaryBatch ~/data/hg19/phyloP46wayAll.bw temp1 3 >tempm
#awk -v name=$2 'BEGIN {printf name}{a=$1+0.1;b=$2+0.1;c=$3+0.1;printf "\t"($1-$3)}END{printf "\n"}' tempm

#paste tempp tempm |awk -v name=$name'_L' 'BEGIN {printf name}{mid=NF/2;for (i=1;i<=5;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (a>0||b>0){printf "\t"; printf log(now);}}}END{printf "\n"}' >temp.dat
#paste tempp tempm |awk -v name=$name'_R' 'BEGIN {printf name}{mid=NF/2;for (i=mid-4;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (a>0||b>0){printf "\t"; printf log(now);}}}END{printf "\n"}' >>temp.dat

#paste tempp tempm |awk 'BEGIN {printf name}{mid=NF/2;for (i=1;i<=5;i++){a=$i;b=$(i+mid);if (a==b) next;now=(a+0.1)/(b+0.1);if (a>0&&b>0){printf "\t"; printf log(now);}}}END{printf "\n"}'

paste tempp tempm |awk 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (a>0||b>0){print i-mid/2,log(now);}}}' >temp.dat
#paste tempp tempm |awk -v name=$2_R 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (a>0||b>0){print i-mid/2,a-b;}}}'
#paste tempp tempm |awk 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (i>1) printf "\t"; printf log(now)} printf "\n"}' >temp.tab

R --no-save --slave --args temp.dat <~/bin/plotXY.R
#R --no-save --slave --args temp.dat <~/bin/plotViolin.R

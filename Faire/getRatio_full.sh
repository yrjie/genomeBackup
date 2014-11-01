if [ $# -lt 2 ]
then
    echo 'Usage: bedfile headNum [bin]'
    exit
fi

bin=20
cell=Gm12878
if [ $# -gt 2 ]
then
    bin=$3
fi

#cell=H1hesc

#pfile=~/Faire/2bp/wgEncodeOpenChromFaire$cell"AlnRep1.plus.bw"
#mfile=~/Faire/2bp/wgEncodeOpenChromFaire$cell"AlnRep1.minus.bw"
#pfile=~/Faire/2bp/AtacGm12878.plus.bw
#mfile=~/Faire/2bp/AtacGm12878.minus.bw
pfile=./FaireGm12878_full.plus.bigWig
mfile=./FaireGm12878_full.minus.bigWig
#pfile=~/Faire/ATAC/AtacGm12878_full.plus.bigWig
#mfile=~/Faire/ATAC/AtacGm12878_full.minus.bigWig

#awk -v win=$win 'BEGIN{OFS="\t"}{print $1,$2-win,$3+win}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=int(($2+$3)/2);print $1,mid-win,mid+win}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{a=$2-win;b=$2+win;print $1,a,b}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=($2+$3)/2;print $1,$2-win,$2}' $1 >temp1
#awk -v win=$win 'BEGIN{OFS="\t"}{mid=($2+$3)/2;print $1,$3,$3+win}' $1 >temp2
#bigWigSummaryBatch ~/Faire/2bp/AtacGm12878.plus.bw temp1 20 >tempp
#bigWigSummaryBatch ~/Faire/2bp/AtacGm12878.minus.bw temp1 20 >tempm

#bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.plus.bw temp1 20 >tempp
#bigWigSummaryBatch ~/Faire/2bp/wgEncodeOpenChromFaireGm12878AlnRep1.minus.bw temp1 20 >tempm

#sort -k7gr $1 |head -n $2 >temp1
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$3-$2}' $1 |sort -k7gr |head -n $2 >temp1
bigWigSummaryBatch $pfile temp1 $bin >tempp
bigWigSummaryBatch $mfile temp1 $bin >tempm


#awk -v win=$win '{if ($3-$2>win) print $0}' $1|sort -k7gr |head -n $2 >temp2
#awk -v win=$win 'BEGIN{OFS="\t"}{a=$2-win/2;b=$2+win/2;print $1,a,b}' temp2 >temp1
##awk -v win=$win 'BEGIN{OFS="\t"}{a=$2;b=$2+win/2;print $1,a,b}' temp2 >temp1
#bigWigSummaryBatch $pfile temp1 $bin >tempp1
#bigWigSummaryBatch $mfile temp1 $bin >tempm1
#
#awk -v win=$win 'BEGIN{OFS="\t"}{a=$3-win/2;b=$3+win/2;print $1,a,b}' temp2 >temp1
##awk -v win=$win 'BEGIN{OFS="\t"}{a=$3-win/2;b=$3;print $1,a,b}' temp2 >temp1
#bigWigSummaryBatch $pfile temp1 $bin >tempp2
#bigWigSummaryBatch $mfile temp1 $bin >tempm2
#
#paste tempp1 tempp2 >tempp
#paste tempm1 tempm2 >tempm

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

paste tempp tempm |awk 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.01)/(b+0.01);if (i>1) printf "\t"; printf log(now)}printf "\n"}' >temp.dat
#paste tempp tempm |awk -v name=$2_R 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (a>0||b>0){print i-mid/2,a-b;}}}'
#paste tempp tempm |awk 'BEGIN{OFS="\t"}{mid=NF/2;for (i=1;i<=mid;i++){a=$i;b=$(i+mid);now=(a+0.1)/(b+0.1);if (i>1) printf "\t"; printf log(now)} printf "\n"}' >temp.tab

#awk -v win=$win '{if (NR==1){for (i=0;i<NF;i++) all[i]=0} for (i=0;i<NF;i++)all[i]+=$(i+1)}END{bin=NF/2; for (i=0;i<bin;i++){printf (i+1)*win/2/bin"\t"all[i]/NR"\t"(i-bin)*win/2/bin"\t"all[i+bin]/NR"\n"}}' temp.dat
awk -v win=$win '{if (NR==1){for (i=0;i<NF;i++) all[i]=0} for (i=0;i<NF;i++)all[i]+=$(i+1)}END{bin=NF/2; for (i=0;i<bin;i++){printf (i+1)*win/bin-win/2"\t"all[i]/NR"\t"(i-bin)*win/bin+win/2"\t"all[i+bin]/NR"\n"}}' temp.dat >tempGP.dat

#echo "call \"plotMean.gp\"" |gnuplot

R --no-save --slave --args temp.dat <plotHeat.R
#R --no-save --slave --args temp.dat <~/bin/plotXY.R
#R --no-save --slave --args temp.dat <~/bin/plotViolin.R

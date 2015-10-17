if [ $# -lt 1 ]
then
    echo 'Usage: alnFile'
    exit
fi

sortFile=sort$RANDOM

LC_ALL=C sort -k4 $1 >$sortFile

awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t.\t"len"\t"preS"\n"$1"\t"$2"\t"$3"\t.\t"(-len)"\t"$6"\n"}}}' $sortFile >aligned.bed

rm $sortFile
#awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t"$1"\t"$2"\t"$3"\t"len"\t"1"\t"1"\n"}}}' top1MRep3.bed >aligned.bedpe
#windowBed -a aligned.bed -b ~/data/hg19/hg19Tss.bed -u -w 0 |awk 'BEGIN{printf "lenAtTSS"}{printf "\t"$4}END{printf "\n"}' >temp.dat
#R --no-save --slave --args temp.dat <~/bin/plotViolin.R
#windowBed -a aligned.bed -b ~/data/hg19/hg19Tss.bed -w 0 |awk 'BEGIN{OFS="\t"}{mid=($2+$3)/2;print mid-$6,$4}' >temp.dat
#R --no-save --slave --args temp.dat <~/bin/plotXY.R
#windowBed -a aligned.bed -b hg19Tss.plus.bed -w 0 |cut -f4 >temp.dat
#windowBed -a aligned.bed -b hg19Tss.plus.bed -w 0 |awk '{print $4+($2+$3)/2-$6}' >temp.dat
#R --no-save --slave --args temp.dat <~/bin/plotHist.R
#awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t"$1"\t"$2"\t"$3"\t"len"\t"1"\t"1"\n"}}}' AtacGm12878AlnRep3_sorted.bed >alignedRep3.bedpe &
#awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3){tmp=preSt;preSt=$2;$2=tmp;tmp=preEd;preEd=$3;$3=tmp}; printf prechr"\t"preSt"\t"preEd"\t"$1"\t"$2"\t"$3"\t"len"\t"1"\t"1"\n"}}}' AtacGm12878AlnRep3_sorted.bed >alignedRep3_swap.bedpe &

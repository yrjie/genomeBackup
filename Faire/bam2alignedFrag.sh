if [ $# -lt 1 ]
then
    echo 'Usage: bamfile'
    exit
fi

tmpBed=tmp$RANDOM

bamToBed -i $1 | LC_ALL=C sort -k4 >$tmpBed
#awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t"len"\t"preS"\n"$1"\t"$2"\t"$3"\t"(-len)"\t"$6"\n"}}}' $tmpBed |sortBed >aligned.bed
grep -v chrM $tmpBed |awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t"len"\t"preS"\n"$1"\t"$2"\t"$3"\t"(-len)"\t"$6"\n"}}}' |sortBed >aligned.bed

rm $tmpBed

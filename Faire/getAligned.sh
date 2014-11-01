if [ $# -lt 1 ]
then
    echo 'Usage: tagBed'
    echo 'the bed file should be sorted by ID'
    exit
fi

awk 'BEGIN{pre="pre";preS="";prechr="";preSt=0;preEd=0;head=end=0;}{if ($6=="+") head=$2; else end=$3; if (!match($4,pre)){pre=substr($4,1,length($4)-1);preS=$6;prechr=$1;preSt=$2;preEd=$3}else {if ($1==prechr&&$6!=preS&&end>head){len=end-head; if (preEd>$3) len=-len; printf prechr"\t"preSt"\t"preEd"\t.\t"len"\t"preS"\n"$1"\t"$2"\t"$3"\t.\t"(-len)"\t"$6"\n"}}}' $1

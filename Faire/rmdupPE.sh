if [ $# -lt 1 ]
then
    echo "Usage: swap.bedpe"
    exit
fi

tmpF=temp$RANDOM.bedpe
#sortBed -i $1 >$tmpF
#awk 'BEGIN{preC="";preS=0;preE=0}{if (preC==$1&&preS==$2&&preE==$6) next; else{preC=$1;preS=$2;preE=$6; print $0}}' $tmpF >temp_rmdup.bedpe
awk 'BEGIN{preC="";preS=0;preE=0}{if (preC==$1&&preS==$2&&preE==$6) next; else{preC=$1;preS=$2;preE=$6; print $0}}' $1 >temp_rmdup.bedpe
rm $tmpF

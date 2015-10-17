if [ $# -lt 1 ]
then
    echo 'Usage: tssSigFile prefix'
    exit
fi

#file=Gm12878TssPro_1000_sig.bed
file=$1
sorted=sorted$RANDOM
prefix=$2

#sort -k9gr $file |centerBed >$sorted
sort -k5gr $file |centerBed >$sorted
#awk 'BEGIN {OFS="\t"}{$9=$7-$8;print $0}' $file |sort -k9gr |centerBed >$sorted
#awk 'BEGIN {OFS="\t"}{$9=$7+$8;print $0}' $file |sort -k9gr |centerBed >$sorted

num=`wc -l $sorted|awk '{print $1}'`
one=`echo $num/10|bc`

split -l $one $sorted splited
suffix=(a b c d e f g h i j)

for i in ${suffix[*]}
do
    echo spliteda$i
#    sh getFragDistr.sh spliteda$i equalTssPro/$prefix$i"p".png
    grep + spliteda$i >spliteda$i"p"
    grep "-" spliteda$i >spliteda$i"m"
    sh getFragDistr.sh spliteda$i"p" equalTssPro/$prefix$i"p".png
    sh getFragDistr.sh spliteda$i"m" equalTssPro/$prefix$i"m".png
    #mv temp.png equalTssPro/$prefix$i.png
done

rm $sorted spliteda*

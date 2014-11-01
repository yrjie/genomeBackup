if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

dat=temp.dat
tmp=tmp$RANDOM

#awk '{print "wget \"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc="$1"\" -q -O html/"$2".html"}' smpLst.txt >getProced.sh

#ls html/*.html |xargs -n1 sh -c 'file=`basename $0`;sh decodeHtml.sh $0 >idVal/${file%%.*}.txt'

#wget "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL7723&id=3825&db=GeoDb_blob41" -O geneAnnot.html

#grep terms= geneAnnot.html >geneAnnot.txt

#python id2symb.py idVal/BRCA1_001.txt geneAnnot.txt geneVal/geneSym.txt
#ls idVal/*.txt|xargs -n1 sh -c 'file=`basename $0`; python id2symb.py $0 geneAnnot.txt >geneVal/gene_$file'


expFile=exp.txt
ls idVal/*.txt |awk -F "." 'NR>1 {printf "\t"}{gsub("idVal/","",$0);printf $1}END{printf "\n"}' >$dat
paste geneVal/geneSym.txt geneVal/gene_*.txt >$expFile

lst=(miR-31 miR-125b-1 miR-125b-2 miR-125a)

for x in ${lst[*]}
do
#    grep \"$x\" $expFile |cut -f3,4,6-|awk 'BEGIN{OFS="\t"}{gsub("\"","",$0);$1=$2"."$1"."NR;print $0}'|cut -f1,3- >>$dat
#    grep \"$x\" $expFile |cut -f3,4,6-|awk 'BEGIN{OFS="\t"}{gsub("\"","",$0);$1=$2"."$1"."NR;a=b=0;for (i=9;i<=14;i++) a+=$i; for (i=15;i<=20;i++) b+=$i; $(NF+1)=a-b; print $0}' |sort -k39gr >>$tmp
    #echo $x
    #grep \"$x\" $expFile |wc -l
    grep "$x" $expFile >>$dat
done
#grep IGF2R $expFile |cut -f4,6- |awk 'BEGIN{OFS="\t"}{gsub("\"","",$1);$1=$1"."NR;print $0}'>>$dat
#grep NESPAS $expFile |cut -f4,6- |awk 'BEGIN{OFS="\t"}{gsub("\"","",$1);$1=$1"."NR;print $0}'>>$dat
#grep AIRN $expFile |cut -f4,6- |awk 'BEGIN{OFS="\t"}{gsub("\"","",$1);$1=$1"."NR;print $0}'>>$dat
#grep SOX9 $expFile |cut -f4,6- |awk 'BEGIN{OFS="\t"}{gsub("\"","",$1);$1=$1"."NR;print $0}'>>$dat

R --no-save --slave --args $dat <plotHeat.R

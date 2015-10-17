if [ $# -lt 1 ]
then
    echo 'Usage: prefix'
    exit
fi

suffix=(abc bc ab cUnique ac bUnique aUnique)

cat $1_abc.bed $1_bc.bed $1_ab.bed $1_cUnique.bed $1_ac.bed $1_bUnique.bed $1_aUnique.bed

#for x in ${suffix[*]}
#do
#    windowBed -a ../pilePeakall3.bed -b $1_$x.bed -w 0 |sort -k4gr|cut -f5- >$1_tmp_$x.bed
#done
#
#cat $1_tmp_abc.bed $1_tmp_bc.bed $1_tmp_ab.bed $1_tmp_cUnique.bed $1_tmp_ac.bed $1_tmp_bUnique.bed $1_tmp_aUnique.bed
#
#rm $1_tmp_*

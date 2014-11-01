if [ $# -lt 1 ]
then
    echo 'Usage: prefix'
    exit
fi

pre=$1
dir=~/chiapet-pipeline-r261/work/$pre

echo '1. generate cpet file'
cut -f1-3,5 $dir/$pre.cpet >$dir/$pre"_cpet.bed"

echo -e "$dir/$pre"_cpet.bed"\tbed4\thg19\t$pre\t$pre"_cpet"\t$pre"_cpet"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >$dir/$pre"_upload.cfg"


echo '2. generate intra-chromosome ipet file'
awk 'BEGIN{OFS="\t"}{gsub(":","\t",$0);gsub("-","\t",$0);} $1==$6&&$11>1 {print $1,$2,$3,$6,$7,$8,$11}' $dir/$pre.cluster >$dir/$pre"_intra.bedpe"

echo -e "$dir/$pre"_intra.bedpe"\tpair.cluster\thg19\t$pre\t$pre"_intra"\t$pre"_intra"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >>$dir/$pre"_upload.cfg"


echo '3. generate inter-chromosome ipet file'
awk 'BEGIN{OFS="\t"}{gsub(":","\t",$0);gsub("-","\t",$0);} $1!=$6&&$11>1 {print $1,$2,$3; print $6,$7,$8}' $dir/$pre.cluster >$dir/$pre"_inter.bed"

echo -e "$dir/$pre"_inter.bed"\tbed3\thg19\t$pre\t$pre"_inter"\t$pre"_inter"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >>$dir/$pre"_upload.cfg"

echo '4. generate cpet_bigWig file'
bedtools genomecov -i $dir/$pre"_cpet.bed" -g ~/genome/human.hg19.genome -bg >$dir/'temp.bedGraph'
bedGraphToBigWig $dir/'temp.bedGraph' ~/genome/human.hg19.genome $dir/$pre"_cpet.bigWig"

echo -e "$dir/$pre"_cpet.bigWig"\tbigWig\thg19\t$pre\t$pre"_cpet_sig"\t$pre"_cpet_sig"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >>$dir/$pre"_upload.cfg"


echo '5. generate cpet peak file'
awk 'BEGIN{OFS="\t"}{gsub(":","\t",$0);gsub("-","\t",$0);print}' $dir/$pre.peak >$dir/$pre"_peak.bed"

echo -e "$dir/$pre"_peak.bed"\tbed4\thg19\t$pre\t$pre"_cpet_peak"\t$pre"_cpet_peak"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >>$dir/$pre"_upload.cfg"


echo '6. generate ipet bigWig file'
awk 'BEGIN{OFS="\t"}{print $1,$2,$2+50;print $5,$6,$6+50}' $dir/$pre.ipet |sortBed >$dir/$pre"_ipet.bed"
bedtools genomecov -i $dir/$pre"_ipet.bed" -g ~/genome/human.hg19.genome -bg >$dir/'temp.bedGraph'
bedGraphToBigWig $dir/'temp.bedGraph' ~/genome/human.hg19.genome $dir/$pre"_ipet.bigWig"

echo -e "$dir/$pre"_ipet.bigWig"\tbigWig\thg19\t$pre\t$pre"_ipet_sig"\t$pre"_ipet_sig"\t{\"assembly\":\"hg19\"}\thttp://biogpu.ddns.comp.nus.edu.sg:8003/" >>$dir/$pre"_upload.cfg"


python ~/bin/uploadTrack.py $dir/$pre"_upload.cfg"


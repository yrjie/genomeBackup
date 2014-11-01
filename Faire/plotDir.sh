if [ $# -lt 1 ]
then
    echo 'Usage: sigfile'
    exit
fi
prefix=${1%%.*}
paste temp $1 |awk 'BEGIN{OFS="\t"}{if (NR%80==0){if ($1=="+") print $0; else{printf $1;for (i=NF;i>=2;i--) printf "\t"$i; printf "\n"}}}' >$prefix"_dir.dat"
~/bin/matrix2png -z -data $prefix"_dir.dat" -mincolor blue -maxcolor green -missingcolor grey -size 4:1  >$prefix"_dir.png"

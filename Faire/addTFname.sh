function writeFile()
{
    awk -v name=$2 'BEGIN{OFS="\t"}{$4=name;print $0}' ~/DNase/data/$1.narrowPeak
}

export -f writeFile
grep Gm12878 ~/DNase/tfbsUnipk.lst |awk 'BEGIN{FS=OFS="\t"}{gsub(" ","_",$2); print $1,$2}' |xargs -n2 sh -c 'writeFile $0 $1'

if [ $# -lt 2 ]
then
    echo 'Usage: nsomePeak TFpeak'
    exit
fi

export nsome=$1
export TF=$2
#ls equalSplit/pilePeaksplit*.bed |xargs -n1 sh -c 'a=`windowBed -a $0 -b $nsome -u -w 200|wc -l`; b=`windowBed -a $0 -b $nsome -u -w 200|windowBed -a stdin -b $TF -u -w 0|wc -l`; c=`sort -k4gr $0|head -n $a|windowBed -a stdin -b $TF -u -w 0|wc -l`; echo -e "$a\t$b\t$c"'
ls pilePeakall*.bed |xargs -n1 sh -c 'echo $0; a=`windowBed -a $0 -b $nsome -u -w 200|wc -l`; b=`windowBed -a $0 -b $nsome -u -w 200|windowBed -b stdin -a $TF -u -w 0|wc -l`; c=`sort -k4gr $0|head -n $a|windowBed -b stdin -a $TF -u -w 0|wc -l`; echo -e "$a\t$b\t$c"'

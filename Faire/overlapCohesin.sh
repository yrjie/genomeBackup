coh=(wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk wgEncodeAwgTfbsHaibGm12878Rad21V0416101UniPk wgEncodeAwgTfbsSydhGm12878Smc3ab9263IggmusUniPk wgEncodeAwgTfbsSydhGm12878Znf143166181apUniPk)
for x in ${coh[@]}
do
    file=~/DNase/data/$x.narrowPeak
    echo $x
    windowBed -a $file -b ~/DNase/data/wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk.narrowPeak -u -w 200|wc -l
    windowBed -a $file -b ~/DNase/data/wgEncodeAwgTfbsSydhGm12878Pol2IggmusUniPk.narrowPeak -u -w 200|wc -l
    windowBed -a $file -b ctcfXpol2.bed -u -w 200|wc -l
done

#$1 TFBS file
#$2 prefix
#$3 cell line
windowBed -a $1 -b ../DNase/Faire/wgEncodeOpenChromFaire$3Pk.narrowPeak -u -w 0 >TFgroup/$2"F.bed"
windowBed -a $1 -b ../DNase/Faire/wgEncodeOpenChromFaire$3Pk.narrowPeak -v -w 5000 >TFgroup/$2"NF.bed"
#windowBed -a $1 -b ../DNase/Faire/wgEncodeOpenChromFaire$3Est10nm30mPk.narrowPeak -u -w 0 >TFgroup/$2"F.bed"
#windowBed -a $1 -b ../DNase/Faire/wgEncodeOpenChromFaire$3Est10nm30mPk.narrowPeak -v -w 5000 >TFgroup/$2"NF.bed"

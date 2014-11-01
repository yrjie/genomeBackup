if [ $# -lt 1 ]
then
    echo 'autorun script'
    exit
fi

ls MACS/*_peaks.bed |xargs -n1 sh -c 'file=`basename $0`;echo $file; windowBed -a ../ATACCd4+/aligned_noM.bed -b $0 -u -w 0 >AtacTag_${file%%_*}.bed'

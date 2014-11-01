if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

windowBed -a /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak -b /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak -v -w 200 |windowBed -a pileall3Filt100_noM_sig.bed -b stdin -u -w 0 >Gm12878Pro.bed

windowBed -b /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak -a /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak -v -w 200 |windowBed -a pileall3Filt100_noM_sig.bed -b stdin -u -w 0 >Gm12878Enh.bed

cat /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak /home/ruijie/data/hg19/broadPeak/wgEncodeBroadHistoneGm12878H3k04me1StdPkV2.broadPeak |windowBed -a pileall3Filt100_noM_sig.bed -b stdin -v -w 0 >Gm12878Other.bed



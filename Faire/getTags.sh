windowBed -a tailPeaks/alltags85-.bed -b tailPeaks/tail60.bed -u -w 0 > tailPeaks/tags60_85-.bed

windowBed -a tailPeaks/alltags85-164.bed -b tailPeaks/tail60.bed -u -w 0 > tailPeaks/tags60_85-164.bed

windowBed -a tailPeaks/alltags164+.bed -b tailPeaks/tail60.bed -u -w 0 > tailPeaks/tags60_164+.bed


windowBed -a tailPeaks/alltags85-.bed -b tailPeaks/tail90.bed -u -w 0 > tailPeaks/tags90_85-.bed

windowBed -a tailPeaks/alltags85-164.bed -b tailPeaks/tail90.bed -u -w 0 > tailPeaks/tags90_85-164.bed

windowBed -a tailPeaks/alltags164+.bed -b tailPeaks/tail90.bed -u -w 0 > tailPeaks/tags90_164+.bed

ls tailPeaks/tags60_* |xargs -n1 sh -c 'echo $0; windowBed -a ../allGm12878TFrmdup.bed -b $0 -u -w 0|wc -l'

ls tailPeaks/tags60_* |xargs -n1 sh -c 'echo $0; windowBed -b ../allGm12878TFrmdup.bed -a $0 -u -w 0|wc -l'


sort -k4g pilePeakall3.bed | head -n600000 >temp.bed

windowBed -a temp.bed -b pilePeakall3Filt100.bed -v -w 0|windowBed -a ../allGm12878TFrmdup.bed -b stdin -u -w 0|windowBed -a stdin -b pilePeakall3Filt100.bed -v -w 0|wc -l

windowBed -a temp.bed -b pilePeakall3Filt100.bed -v -w 0|windowBed -a ../allGm12878TFrmdup.bed -b stdin -u -w 0|windowBed -a stdin -b pilePeakall3Filt100.bed -v -w 0|windowBed -a temp.bed -b stdin -u -w 0  >temp1.bed

windowBed -a temp.bed -b pilePeakall3Filt100.bed -v -w 0|windowBed -a stdin -b ../allGm12878TFrmdup.bed -v -w 0 >temp2.bed


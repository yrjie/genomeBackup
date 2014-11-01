sh getRocPeakShuff.sh fragLen/pilePeakall3.bed Gm12878 allGm12878TFrmdup_merge.bed pileall3
sh getRocPeakShuff.sh openEM/discrete2dAllth_tfBin3.bed Gm12878 allGm12878TFrmdup_merge.bed discrete2dAllth_tfBin3
sh getRocPeakShuff.sh openEM/discrete2dAllth_tfBin5.bed Gm12878 allGm12878TFrmdup_merge.bed discrete2dAllth_tfBin5
sh getRocPeakShuff.sh openEM/discrete2dAllth_tfBin8.bed Gm12878 allGm12878TFrmdup_merge.bed discrete2dAllth_tfBin8
sh getRocPeakShuff.sh openEM/discrete2dAllth_tfBin15.bed Gm12878 allGm12878TFrmdup_merge.bed discrete2dAllth_tfBin15

cut -f1-3,5 MACS/AtacGm12878_peaks.bed >temp1.bed
sh getRocPeakShuff.sh temp1.bed Gm12878 allGm12878TFrmdup_merge.bed macs

#windowBed -a AtacGm12878_lowSig.bed -b allCtcf.bed -u -w 0 |grep -P "chr1\t" |head -n 5000 >pos.bed
windowBed -a AtacGm12878_lowSig.bed -b uncovTF.bed -u -w 0 |grep -P "chr1\t" |head -n 5000 >pos.bed
#windowBed -a AtacGm12878_lowSig.bed -b ../allGm12878TF.bed -v -w 200 |grep -P "chr1\t" |head -n 5000 >neg.bed
windowBed -a AtacGm12878_lowSig.bed -b ../allGm12878TF.bed -v -w 200 |windowBed -a stdin -b ../allGm12878Histone.bed -v -w 200|grep -P "chr1\t" |head -n 5000 >neg.bed


ls AtacPk_*.bed | xargs -n1 sh -c 'windowBed -a ../../ATAC/alignedRep3_sorted_filt100.bed -b $0 -u -w 0 > AtacTag${0#AtacPk*}'

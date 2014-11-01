ell3=peak/ell3Raw_filtered.bed
ERVK=mm9ERVK.bed
TF=$1
wc -l $ell3
wc -l $ERVK
wc -l $TF
windowBed -a $ell3 -b $ERVK -u -w 0 |wc -l
windowBed -a $ell3 -b $TF -u -w 500 |wc -l
windowBed -a $TF -b $ERVK -u -w 0 |wc -l
windowBed -a $TF -b $ERVK -u -w 0 |windowBed -a stdin -b repeat/ell3Raw_filtered_ERVK.bed -u -w 500|wc -l

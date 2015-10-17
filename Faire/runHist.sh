cut -f1-4 AtacPkTag_pileall3_noM.bed |uniq -c |awk '{print $1}' >temp.dat

R --no-save --slave --args temp.dat < plotTagHist.R



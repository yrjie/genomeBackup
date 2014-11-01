paste Gm12878Tss.bed templ |awk '{if ($9>-15&&$9<0&&$10>30) print $0}' >temp0
paste Gm12878Tss.bed templ |awk '{if ($9>0&&$9<15&&$10>30) print $0}' >temp2
paste Gm12878Tss.bed templ |awk '{if (($9>15||$9<-15)&&$10>30) print $0}' >temp1
wc -l temp0
wc -l temp1
wc -l temp2
sh getDiffRNA.sh temp0 temp1 temp2

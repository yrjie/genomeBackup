#ls equalSplit/pilePeaksplit*_3.bed |xargs -n1 sh -c 'echo $0; cut -f4 $0 >temp.dat; R --no-save --slave --args temp.dat <~/bin/plotHist.R'

ls equalSplit/pilePeaksplit*.bed |xargs -n1 sh -c 'echo $0;python checkTFRatio.py $0 ../allGm12878TF.bed |grep CTCF_w|cut -f2'

#ls equalSplit/pilePeaksplit*_0.bed|xargs -n1 awk 'BEGIN{sum=0}{sum+=($3-$2)*$4}END{print sum}'

#ls equalSplit/pilePeaksplit*_0.bed |xargs -n1 sh -c 'echo $0; sh ~/bin/getLenDistr.sh $0'

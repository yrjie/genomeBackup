if [ $# -lt 1 ]
then
    echo 'Usage: datfile'
    exit
fi

a=("normal" "gamma" "poisson" "\"negative binomial\"" "exponential")

#for ((i=0;i<${#a[*]};i++))
#do
#    echo ${a[$i]}
#    R --no-save --slave --args $1 ${a[$i]} <~/bin/plotFit.R
#done

R --no-save --slave --args $1 "normal" <~/bin/plotFit.R
R --no-save --slave --args $1 "gamma" <~/bin/plotFit.R
R --no-save --slave --args $1 "poisson" <~/bin/plotFit.R
R --no-save --slave --args $1 "negative binomial" <~/bin/plotFit.R
R --no-save --slave --args $1 "exponential" <~/bin/plotFit.R
#R --no-save --slave --args $1 "beta" <~/bin/plotFit.R

#awk '$4>=0&&$4<0.090909091 {print $0}' filt100_bratio.bed | windowBed -a stdin -b AtacPkTagFilt100_pileall3_noM.bed -w 0|cut -f1-3|uniq -c |awk '{print $1}' >temp1.dat

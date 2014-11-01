if [ $# -lt 3 ]
then
    echo 'Usage: pos.bed neg.bed align.bed'
    exit
fi

#windowBed -a $3 -b $1 -u -w 0 |awk 'BEGIN{a=b=c=d=0;OFS="\t"}{if ($4<0) x=-$4; else x=$4; if (x<150) a++; else if (x<250) b++; else if(x<400) c++; else d++;}END{print a/NR,b/NR,c/NR,d/NR}' >pos.tab
#windowBed -a $3 -b $2 -u -w 0 |awk 'BEGIN{a=b=c=d=0;OFS="\t"}{if ($4<0) x=-$4; else x=$4; if (x<150) a++; else if (x<250) b++; else if(x<400) c++; else d++;}END{print a/NR,b/NR,c/NR,d/NR}' >neg.tab
#tab2arff.py pos.tab neg.tab ATAC >temp.arff

#windowBed -a temp.bed -b uncovCtcf.bed -u -w 0|awk 'BEGIN{OFS="\t"}{print $1,$2+$10-50,$2+$10+50,$4}' >pos.bed
#windowBed -a temp.bed -b uncovCtcf.bed -v -w 200|awk 'BEGIN{OFS="\t"}{print $1,$2+$10-50,$2+$10+50,$4}' >neg.bed

windowBed -a $1 -b $3 -w 0 >overlap.tmp
#centerBed $1 |slopBed -g ~/genome/human.hg19.genome -b 50| windowBed -a stdin -b $3 -w 0 >overlap.tmp
#python getLenTab.py overlap.tmp >pos.tab
#python getLMR.py overlap.tmp >pos.tab
python getFuzz.py overlap.tmp >pos.tab
#R --no-save --slave --args pos.tab <~/bin/plotMeanCov.R
windowBed -a $2 -b $3 -w 0 >overlap.tmp
#centerBed $2 |slopBed -g ~/genome/human.hg19.genome -b 50| windowBed -a stdin -b $3 -w 0 >overlap.tmp
#python getLenTab.py overlap.tmp >neg.tab
#python getLMR.py overlap.tmp >neg.tab
python getFuzz.py overlap.tmp >neg.tab
#R --no-save --slave --args neg.tab <~/bin/plotMeanCov.R
tab2arff.py pos.tab neg.tab ATAC >temp.arff

#sh tab2Box.sh pos.tab 150_pos
#sh tab2Box.sh neg.tab 150_neg

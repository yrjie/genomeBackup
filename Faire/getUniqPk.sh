if [ $# -lt 5 ]
then
    echo 'Usage: 5 peakfiles'
    exit
fi

windowBed -a $1 -b $2 -v -w 0 |windowBed -a stdin -b $3 -v -w 0|windowBed -a stdin -b $4 -v -w 0|windowBed -a stdin -b $5 -v -w 0 >temp1.bed

windowBed -a $2 -b $1 -v -w 0 |windowBed -a stdin -b $3 -v -w 0|windowBed -a stdin -b $4 -v -w 0|windowBed -a stdin -b $5 -v -w 0 >temp2.bed

windowBed -a $3 -b $1 -v -w 0 |windowBed -a stdin -b $2 -v -w 0|windowBed -a stdin -b $4 -v -w 0|windowBed -a stdin -b $5 -v -w 0 >temp3.bed

windowBed -a $4 -b $1 -v -w 0 |windowBed -a stdin -b $2 -v -w 0|windowBed -a stdin -b $3 -v -w 0|windowBed -a stdin -b $5 -v -w 0 >temp4.bed

windowBed -a $5 -b $1 -v -w 0 |windowBed -a stdin -b $2 -v -w 0|windowBed -a stdin -b $3 -v -w 0|windowBed -a stdin -b $4 -v -w 0 >temp5.bed



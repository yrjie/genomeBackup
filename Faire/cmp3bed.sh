err(){
    echo "[`date +'%Y-%m-%dT%H:%M:%S%z'`]: $*" >&2
}

if [ $# -lt 3 ]
then
#    echo 'Usage: bed1 bed2 targetBed'
    err 'Usage: bed1 bed2 targetBed'
    exit
fi

wc -l $1
wc -l $2
wc -l $3

echo "1 vs 2"
windowBed -a $1 -b $2 -u -w 0|wc -l
windowBed -b $1 -a $2 -u -w 0|wc -l

echo "1 vs 3"
windowBed -a $1 -b $3 -u -w 0|wc -l
windowBed -b $1 -a $3 -u -w 0|wc -l

echo "2 vs 3"
windowBed -a $2 -b $3 -u -w 0|wc -l
windowBed -b $2 -a $3 -u -w 0|wc -l

echo "1 vs 2 vs 3"
windowBed -a $1 -b $2 -u -w 0|windowBed -a stdin -b $3 -u -w 0|wc -l
windowBed -a $1 -b $2 -u -w 0|windowBed -b stdin -a $3 -u -w 0|wc -l
windowBed -b $1 -a $2 -u -w 0|windowBed -a stdin -b $3 -u -w 0|wc -l
windowBed -b $1 -a $2 -u -w 0|windowBed -b stdin -a $3 -u -w 0|wc -l

#echo "(1 no 2) vs 3"
#windowBed -a $1 -b $2 -v -w 0|windowBed -a stdin -b $3 -u -w 0|wc -l
#windowBed -a $1 -b $2 -v -w 0|windowBed -b stdin -a $3 -u -w 0|wc -l
#
#echo "(2 no 1) vs 3"
#windowBed -a $2 -b $1 -v -w 0|windowBed -a stdin -b $3 -u -w 0|wc -l
#windowBed -a $2 -b $1 -v -w 0|windowBed -b stdin -a $3 -u -w 0|wc -l
#
#echo "(1 no 3) vs 2"
#windowBed -a $1 -b $3 -v -w 0|windowBed -a stdin -b $2 -u -w 0|wc -l
#windowBed -a $1 -b $3 -v -w 0|windowBed -b stdin -a $2 -u -w 0|wc -l

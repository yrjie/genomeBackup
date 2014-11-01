if [ $# -lt 2 ]
then
    echo 'Usage: peak1 peak2'
    echo 'The 4th column is the score'
    exit
fi

tmpFile=temp$RANDOM.dat
windowBed -a $1 -b $2 -w 0|cut -f4,8 |shuf -n200000 >$tmpFile
R --no-save --slave --args $tmpFile <~/bin/plotXY.R
rm $tmpFile

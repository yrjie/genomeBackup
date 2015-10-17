if [ $# -lt 2 ]
then
    echo 'Usage: predicted answer'
    exit
fi

grep "+" $1 |wc -l
grep "+" $1 |windowBed -a $2 -b stdin -u -w 0|grep "+" |wc -l

grep "\-" $1 |wc -l
grep "\-" $1 |windowBed -a $2 -b stdin -u -w 0|grep "-" |wc -l

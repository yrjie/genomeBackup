if [ $# -lt 2 ]
then
    echo 'Usage: XYdata binNum'
    exit
fi

file=$RANDOM.dat
python ~/bin/binData.py $1 $2 >$file
R --no-save --slave --args $file <~/bin/plotBox.R
rm $file

echo $1 $2 $3 \'$1\'

echo $# $*
echo $@
if [ $# -lt 2 ]
then
    echo "lt 2"
else
    echo "not lt 2"
fi
ls $2

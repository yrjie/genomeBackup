if [ $# -lt 1 ]
then
    echo 'Usage: # of run'
    exit
fi

for ((i=0;i<=$1;i++)) do
    python EMtest.py
done

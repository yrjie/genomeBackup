if [ $# -lt 1 ]
then
    echo 'Usage: prefixTarget [prefixControl]'
    echo 'output to stdout'
    exit
fi

ls $1*

if [ $# -gt 1 ]
then
    echo ===
    ls $2*
fi


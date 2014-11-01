#!/bin/bash

if [ $# -lt 2 ]
then
    echo "Usage: bigWig bed"
    exit
fi

tmp1=/tmp/temp1$RANDOM
tmp2=/tmp/temp2$RANDOM
tmp3=/tmp/temp3$RANDOM
bigWigSummaryBatch $1 $2 1 >$tmp1
bigWigSummaryBatch $1 $2 1 -type=coverage >$tmp2
awk '{print $3-$2}' $2 >$tmp3
paste $tmp1 $tmp2 $tmp3 |awk '{if ($2>0) print $1*$2*$3; else print 0}'

rm $tmp1 $tmp2 $tmp3

if [ $# -lt 4 ]
then
    echo 'Usage: 85-.bed 85-164.bed 164+.bed TFfile'
    exit
fi

tmpFile=temp$RANDOM.bed

echo a unique
windowBed -a $1 -b $2 -v -w 0|windowBed -a stdin -b $3 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_aUnique.bed

echo b unique
windowBed -a $2 -b $1 -v -w 0|windowBed -a stdin -b $3 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_bUnique.bed

echo c unique
windowBed -a $3 -b $1 -v -w 0|windowBed -a stdin -b $2 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_cUnique.bed

echo ab
windowBed -a $1 -b $2 -u -w 0|windowBed -a stdin -b $3 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_ab.bed

echo ac
windowBed -a $1 -b $3 -u -w 0|windowBed -a stdin -b $2 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_ac.bed

echo bc
windowBed -a $2 -b $3 -u -w 0|windowBed -a stdin -b $1 -v -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_bc.bed

echo abc
windowBed -a $1 -b $2 -u -w 0|windowBed -a stdin -b $3 -u -w 0 >$tmpFile
wc -l $tmpFile
windowBed -a $4 -b $tmpFile -u -w 0|wc -l
windowBed -b $4 -a $tmpFile -u -w 0|wc -l
cp $tmpFile allComb_abc.bed

rm $tmpFile

if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

python openTSS.py AtacPkTag_Gm12878TssProTopAtac.bed temp.bed 1 1 >temp.dat
python ~/bin/runLineReg.py temp.dat 1
awk 'BEGIN{OFS="\t"}{a=8.52;b=-395.44;c=-142.736;print a*$1+b*$2+c*$3,$4}' temp.dat >temp1.dat
R --no-save --slave --args temp1.dat <~/bin/plotXY.R


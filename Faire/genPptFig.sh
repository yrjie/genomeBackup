if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

python openEM.py AtacPkTagFilt100_pileall3_noM.bed 8 temp.bed 1 >openEM_ppt.log
tail -n+2054364 openEM_ppt.log |sort -k1gr |less

head -n684790 openEM_ppt.log >iter1.log

head -n684790 openEM_ppt.log |sort -k1gr |less

head -n684790 openEM_ppt.log |awk 'NF==3&&$3>0.5{print $1}' >temp.dat

head -n684790 openEM_ppt.log |awk 'NF==3&&$3>0.5{print $2}' >temp.dat

head -n684790 openEM_ppt.log |awk 'NF==3&&$3<0.3&&$1>7{print $1}' >temp.dat


awk 'NF==3&&$3>0.5{print $1}' iter2.log >temp.dat

awk 'NF==3&&$3<0.9999&&$1>7{print $2}' iter3.log >temp.dat
awk 'NF==3&&$3>0.9999&&$1>7{print $2}' iter3.log >temp.dat

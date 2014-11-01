if [ $# -lt 1 ]
then
    echo 'Usage: histoneFiles'
    exit
fi

#ls $1 |xargs -n1 sh -c 'a=${0##*Gm12878};a=${a%.*};b=`windowBed -a ../ATAC/uncovTFrmdup_pile.bed -b $0 -u -w 0|wc -l`; echo -e $a"\t"$b' >temp.txt
#ls $1 |xargs -n1 sh -c 'a=${0##*Gm12878};a=${a%.*};b=`windowBed -a ../ATAC/covTFrmdup_pile.bed -b $0 -u -w 0|wc -l`; echo -e $a"\t"$b' >temp.txt
#ls $1 |xargs -n1 sh -c 'a=${0##*Gm12878};a=${a%.*};b=`windowBed -a ../shuffled.bed -b $0 -u -w 0|wc -l`; echo -e $a"\t"$b' >temp.txt

ls $* |xargs -n1 sh -c 'a=${0##*Gm12878};a=${a%.*}; echo $a'

echo uncovTF
ls $* |xargs -n1 sh -c 'windowBed -a ../ATAC/uncovTFrmdup_pile.bed -b $0 -u -w 0|wc -l'
echo covTF
ls $* |xargs -n1 sh -c 'windowBed -a ../ATAC/covTFrmdup_pile.bed -b $0 -u -w 0|wc -l'
echo shuffled
ls $* |xargs -n1 sh -c 'windowBed -a ../shuffled.bed -b $0 -u -w 0|wc -l'

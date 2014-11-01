
#ls test/AtacGm12878*NoStL100Ex* |xargs -n1 sh -c 'echo $0;sort -k7gr $0|head -n100000|windowBed -a stdin -b ../allGm12878TF.bed -u -w 0|wc -l'
ls test/AtacGm12878Rmdup*St.bed |xargs -n1 sh -c 'echo $0;sort -k7gr $0|head -n100000 >temp.bed; windowBed -a temp.bed -b ../allGm12878TF.bed -u -w 0|wc -l;windowBed -a temp.bed -b ../shuffledTF.bed -u -w 0|wc -l;  windowBed -b temp.bed -a ../allGm12878TF.bed -u -w 0|wc -l; windowBed -b temp.bed -a ../shuffledTF.bed -u -w 0|wc -l;'

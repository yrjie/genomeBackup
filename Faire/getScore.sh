windowBed -b FaireHela3_D.bed -a temp -u -w 0 |cut -f5|awk '{if (NR==1) printf "Hela_FD"; printf "\t"$1}END{printf "\n"}' >temp.dat
windowBed -b FaireHela3_noD.bed -a temp -u -w 0 |cut -f5|awk '{if (NR==1) printf "Hela_FnoD"; printf "\t"$1}END{printf "\n"}' >>temp.dat
R --no-save --slave --args temp.dat <~/bin/plotBox.R

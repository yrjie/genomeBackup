awk 'BEGIN{OFS="\t"}{l=0;r=0;for (i=2;i<7;i++)l+=$i;for (i=13;i<18;i++)r+=$i;print l-r,$0}' Gm12878FaireSp1_c1.sig|sort -k1 -g >temp.sig

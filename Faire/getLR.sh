prefix=$1
if [ $2 == "sig" ]
then
    paste $prefix.sig $prefix.bed|awk 'BEGIN{OFS="\t"}{l=0;r=0;for (i=6;i<9;i++)l+=$i;for (i=13;i<16;i++)r+=$i;if(l>r) print $21,$22,$23,l-r}' >templ.bed
    paste $prefix.sig $prefix.bed|awk 'BEGIN{OFS="\t"}{l=0;r=0;for (i=6;i<9;i++)l+=$i;for (i=13;i<16;i++)r+=$i;if(l<r) print $21,$22,$23,l-r}' >tempr.bed
else
    paste $prefix.phyl $prefix.bed|awk 'BEGIN{OFS="\t"}{l=0;r=0;for (i=5;i<10;i++)l+=$i;for (i=12;i<17;i++)r+=$i;if(l>r) print $21,$22,$23,l-r}' >templ.bed
    paste $prefix.phyl $prefix.bed|awk 'BEGIN{OFS="\t"}{l=0;r=0;for (i=5;i<10;i++)l+=$i;for (i=12;i<17;i++)r+=$i;if(l<r) print $21,$22,$23,l-r}' >tempr.bed
fi
wc -l templ.bed
wc -l tempr.bed
closestBed -a templ.bed -b ../../data/knownGene.bed -t first |cut -f10|sort -r|uniq -c
closestBed -a tempr.bed -b ../../data/knownGene.bed -t first |cut -f10|sort -r|uniq -c

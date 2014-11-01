ls ../data/Gm12878Egr1_Fclust*.bed|xargs -n1 -I xx sh -c "echo xx; closestBed -a xx -b ../data/knownGene.bed -t first -d|awk '{if (\$NF<1000) print \$0}' |cut -f12|sort|uniq -c"

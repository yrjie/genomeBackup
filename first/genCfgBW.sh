ls /tmp/*raph|xargs -n1 sh -c 'bedGraphToBigWig $0 ~/genome/human.hg19.genome ${0%%.*}.bigWig'

ls /tmp/*.bigWig|awk 'BEGIN{OFS="\t"}{label=$0;gsub("/tmp/","",label);gsub(".bigWig","",label);print $0,"bigWig", "hg19", "AR-ERG",label,label, "{\"assembly\":\"hg19\"}","http://biogpu.ddns.comp.nus.edu.sg:8003/" }'

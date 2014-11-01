#java -cp config:classes:lib/* -Xmx64m -Dlog4j.configuration=log4j.properties -Djppf.config=jppf.properties -Djava.util.logging.config.file=config/logging.properties org.jppf.application.template.TemplateApplicationRunner

#java -Djava.awt.headless=false -cp "config:bin:lib/*" PeakClassifier -dataDir /home/ruijie/ELL3/Tcf3_clust/ -print -sumNorm -medianNorm -sort -winsize 100000 -clusternum 5 -JPPFmode 2 -peakfile1 $*

java -Dlog4j.configuration=config/log4j.properties -Djava.awt.headless=false -cp "config:bin:lib/*" PeakClassifier -dataDir ~/Faire/H1hesc/ -print -sumNorm -medianNorm -sort -winsize 100000 -clusternum 5 -JPPFmode 2 -peakfile1 ./ctcf_test1.bed

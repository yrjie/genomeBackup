infile='DHMS.region'
#infile='data/LD_H3K4me3_786-O_138_a_summits.bed'
#infile='data/138vsNS_peaks.bed'
#infile='MAnorm_Linux_R_Package/MAnorm_result.xls'
#arrayFile='microarray_folder/JARID_mas5.txt'
arrayFile='microarray_folder/JARID_rma.txt'
python getMAPlot.py $infile microarray_folder/U133P2_gene.xls $arrayFile
#python getMAPlot.py $infile microarray_folder/affyU133P2.xls $arrayFile
R --no-save --slave <rplot.R

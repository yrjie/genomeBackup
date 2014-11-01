args=commandArgs()
if (length(args)<5){
	print("Usage: outfile")
	print("cd to the dir containing CEL files")
	q()
}

#http://homer.salk.edu/homer/basicTutorial/affymetrix.html

#source("http://bioconductor.org/biocLite.R")

# mouse
#biocLite("mogene20sttranscriptcluster.db")

# human
#biocLite("hugene21sttranscriptcluster.db")

outfile=args[5]

library(affy)
#library(mogene20sttranscriptcluster.db)
library(hgu133plus2.db)
data <- ReadAffy()
eset <- rma(data)
#write.exprs(eset,file=outfile)
my_frame <- data.frame(exprs(eset))
#Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))

Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))

all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
write.table(all,file=outfile,sep="\t")

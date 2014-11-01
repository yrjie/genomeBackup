args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [infile2] <plotXY.R")
q()
}
library(affy)
infile=args[5]
outfile='temp.png'
dat=as.matrix(read.table(infile,sep='\t'))
x=dat[,1]
if (length(args)>=6){
	y=as.matrix(read.table(args[6],sep='\t'))[,1]
} else y=dat[,2]

#library(preprocessCore)
#norm=normalize.quantiles(cbind(x,y))
#x=norm[,1]
#y=norm[,2]

png(outfile)
ma.plot(rowMeans(log2(cbind(x,y))),log2(y)-log2(x),cex=1)
title('MA plot')
dev.off()

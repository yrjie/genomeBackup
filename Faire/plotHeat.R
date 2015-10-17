args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [sumNorm] [rowNorm] <plotHeat.R")
q()
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- qnt[1]-H
  y[x > (qnt[2] + H)] <- qnt[2]+H
  y
}

library(gplots)
infile=args[5]
outfile='temp.png'
dat=as.matrix(read.table(infile,header=F))
if ("sumNorm" %in% args)
	dat=sweep(dat,2,colSums(dat),FUN="/")
#print(length(dat[dat<mean(dat)]))
#dat=(dat-mean(dat))/sd(dat)
#dat=remove_outliers(dat)
dDat=dim(dat)
print(dDat)
png(outfile, width=max(1200,30*dDat[2]), height=max(1200,30*dDat[2]))
sc="none"
if ("rowNorm" %in% args)
	sc="row"
#heatmap.2(dat, dendrogram="col", col=bluered(256), scale=sc, key=T, keysize=0.8, density.info="none", trace="none",cexCol=2.3, labRow=NA, margins=c(20,5),Rowv=FALSE)
#heatmap.2(dat, dendrogram="row", col=bluered(256), scale=sc, key=T, keysize=0.8, density.info="none", trace="none",cexCol=2.3, labRow=NA, margins=c(20,5),Colv=FALSE)
heatmap.2(dat, dendrogram="none", col=bluered(256), scale=sc, key=T, keysize=0.8, density.info="none", trace="none",cexCol=2.3, labRow=NA, labCol=NA, margins=c(20,5),Colv=FALSE, Rowv=FALSE)
#clustering: Colv Rowv
dev.off()

Args<-commandArgs()
if (length(Args)<5){
	print("Usage: R --no-save --slave --args infile [outfile]")
	q()
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

outfile='temp.png'
if (length(Args)>5){
	outfile=Args[6]
}

data <- read.table(Args[5], header=FALSE, sep ="\t")
data=as.numeric(as.matrix(data))
print(summary(data))
data=remove_outliers(data)
#png(paste(gsub("\\$","_",Args[5]),".png",sep=""))
png(outfile)
opar=par(ps=18)
hist(data,breaks=50,main="",col='red',xlab="x value", ylab="Counting")
if (length(Args)>5){
	mtext(paste(outfile),side=3,line=0,cex=2)
}
dev.off()

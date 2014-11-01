args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [infile2] <plotXY.R")
q()
}
infile=args[5]
outfile='temp.png'

remove_outliers <- function(x, y, na.rm = TRUE, ...) {
  qntx <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  qnty <- quantile(y, probs=c(.25, .75), na.rm = na.rm, ...)
  Hx <- 10 * IQR(x, na.rm = na.rm)
  Hy <- 10 * IQR(y, na.rm = na.rm)
  xrm <- x
  xrm[x < (qntx[1] - Hx)] <- NA
  xrm[x > (qntx[2] + Hx)] <- NA
  xrm[y < (qnty[1] - Hy)] <- NA
  xrm[y > (qnty[2] + Hy)] <- NA
  xrm=xrm[!is.na(xrm)]
  xrm
}

dat=as.matrix(read.table(infile,sep='\t'))
x=dat[,1]
if (length(args)>=6){
	y=as.matrix(read.table(args[6],sep='\t'))[,1]
} else y=dat[,2]
x=as.numeric(x)
y=as.numeric(y)
xrm=remove_outliers(x, y)
yrm=remove_outliers(y, x)
png(outfile)
#plot(xrm,yrm,xlab='regionVal',ylab='validation',lty=1)
plot(x,y,xlab='regionVal',ylab='validation',lty=1)
title('validation-regionVal correlation')
#corr=cor(xrm,yrm)
corr=cor(x,y)
print(corr)
mtext(paste('correlation=',corr,sep=''),side=3,line=0,cex=2)
dev.off()

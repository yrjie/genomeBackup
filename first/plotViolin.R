args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [name] <plotViolin.R")
q()
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.15, .85), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- qnt[1] - H
  y[x > (qnt[2] + H)] <- qnt[2] + H
  y =y[!is.na(y)]
  y
}

library(vioplot)
infile=args[5]
name=""
if (length(args)>=6){
	name=args[6]
}
outfile='temp.png'
fi=file(infile,open="r")
dat=list()
while (length(line<-readLines(fi,n=1,warn=FALSE))>0){
	vec=strsplit(line,"\t")
	lenv=length(vec[[1]])
	cname=vec[[1]][1]
	vec=list(as.numeric(vec[[1]][2:lenv]))
#	dat[cname]=vec
	print(summary(vec[[1]]))
	dat[cname]=list(remove_outliers(vec[[1]]))
	print(length(dat[[cname]]))
}
nameVector=names(dat)
names(dat)[1]='x'
print(nameVector)
png(outfile)
#boxplot(dat,main=name,range=0)
do.call(vioplot,c(dat,list(names=nameVector,col="cyan")))
dev.off()

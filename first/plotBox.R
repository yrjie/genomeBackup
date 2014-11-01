args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [name] <plotBox.R")
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
}
print(names(dat))
png(outfile, width=max(450,120*length(dat)), height=max(450,120*length(dat)),units="px")
#boxplot(dat,main=name,range=0)
boxplot(dat,main=name)
dev.off()

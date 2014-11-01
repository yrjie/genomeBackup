args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [name] [rowNorm] <plotMeanCov.R")
q()
}

dat=as.matrix(read.table(args[5]))
x0=dim(dat)[2]
x=(-x0/2):(x0/2-1)
name=""
if (length(args)>5){
	name=paste(args[6])
}
png("temp.png")
if ("rowNorm" %in% args)
	dat=t(apply(dat,1,function(x)(if (max(x)-min(x)<1e-8){x-min(x)}else (x-min(x))/(max(x)-min(x)))))
#print(apply(dat,2,mean))
plot(x,apply(dat,2,mean),xlab="pos",ylab="Mean coverage",type='l',lwd=4, main=name)
dev.off()

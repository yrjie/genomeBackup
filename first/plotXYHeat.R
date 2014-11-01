args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [infile2] <plotXY.R")
q()
}
infile=args[5]
outfile='temp.png'
dat=as.matrix(read.table(infile,sep='\t'))
x=dat[,1]
if (length(args)>=6){
	y=as.matrix(read.table(args[6],sep='\t'))[,1]
} else y=dat[,2]
df <- data.frame(x=x,y=y)
png(outfile)
k <- with(df,MASS:::kde2d(x,y))
filled.contour(k,xlab='regionValue',ylab='validation', color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')))
corr=cor(x,y)
print(corr)
mtext(paste('correlation=',corr,sep=''),side=3,line=0,cex=2)
dev.off()

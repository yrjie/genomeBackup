args=commandArgs()
if (length(args)<5){
print("Usage: R --no-save --slave --args infile [name] <plotBox.R")
q()
}

infile=args[5]
name=""
if (length(args)>=6){
	name=args[6]
}
outfile='temp.png'
dat=as.matrix(read.table(infile,sep='\t'))
x=dat[,1]
if (length(args)>=6){
        y=as.matrix(read.table(args[6],sep='\t'))[,1]
} else y=dat[,2]

y=as.double(y)
ylimit = max(y)*1.1

png(outfile, width=max(450,30*length(x)), height=max(450,30*length(x)),units="px")
par(mar=c(15,5,5,2))
bp <- barplot(y,names.arg=x,ylim=c(0,ylimit),cex.names=1.5,col=rainbow(6),las=2)
#bp <- barplot(y,names.arg=x,ylim=c(0,ylimit),cex.names=1.5,col='red',las=2)
dev.off()

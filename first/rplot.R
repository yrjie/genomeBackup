infile='MAPlot.dat'
#outfile='ChipDiffPlot.png'
#outfile='MACSPlot.png'
outfile='MAPlot.png'
dat=as.matrix(read.table(infile,sep='\t'))
x=dat[,1]
y=dat[,2]
png(outfile)
#plot(x,y,xlab='confidence',ylab='log(expression fold)',lty=1)
#title('expression-ChipDiff correlation')
#plot(x,y,xlab='MACS',ylab='log(expression fold)',lty=1)
#title('expression-MACS correlation')
plot(x,y,xlab='regionVal',ylab='log(expression fold)',lty=1)
title('expression-regionVal correlation')
corr=cor(x,y)
print(corr)
mtext(paste('correlation=',corr,sep=''),side=3,line=0)
#y=sin(x)
#plot(x,y,ylab='y1',col=mycol[2],lty=2)
dev.off()

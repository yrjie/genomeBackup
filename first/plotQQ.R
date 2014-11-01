Args<-commandArgs()
if (length(Args)<5){
	print("Usage: R --no-save --slave --args infile")
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

data <- read.table(Args[5], header=FALSE, sep ="\t")
print(summary(data[,1]))
data=remove_outliers(data[,1])
#data=data[,1]
#png(paste(gsub("\\$","_",Args[5]),".png",sep=""))
png('temp.png')
qqnorm(data)
qqline(data)
dev.off()

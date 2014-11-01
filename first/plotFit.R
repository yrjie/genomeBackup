Args<-commandArgs()
if (length(Args)<6){
	print("Usage: R --no-save --slave --args infile type [outfile]")
	print('Type should be one of: "beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial", "normal", "poisson", "t" and "weibull"')
        q()
}

outfile="temp.png"
if (length(Args)>=7)
#	outfile=paste(gsub(".",Args[7],"_"),".png",sep="")
	outfile=Args[7]

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
#  y[x < (qnt[1] - H)] <- qnt[1] - H
#  y[x > (qnt[2] + H)] <- qnt[2] + H
  y=y[x >= (qnt[1] - H) & x <= (qnt[2] + H)]
  y
}

x=as.numeric(as.matrix(read.table(Args[5])))
print(paste("shift",min(x)-1,sep=" "))
x=x-min(x)+1
x=remove_outliers(x)

type=Args[6]
library(MASS)
#x=rgamma(100,shape=3.5,scale=0.5)
#fitdistr(x,"gamma")
#fitdistr(x,densfun=dweibull,start=list(scale=1,shape=2))
#fitdistr(x,"normal")
fit=fitdistr(x,type)
#fit=fitdistr(x,densfun=dchisq, start=list(df=1,ncp=2))
print(fit)

h<-hist(x,breaks=20)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(x),max(x),length=100)
yfit=0
ytest=0
ntest=10000

# poisson and nb are integer-based, not good

if (length(grep(type,"normal"))>0){
	yfit<-dnorm(xfit,mean=fit$estimate[1],sd=fit$estimate[2])
	ytest<-rnorm(ntest,mean=fit$estimate[1],sd=fit$estimate[2])
}
if (length(grep(type,"gamma"))>0){
	yfit<-dgamma(xfit,shape=fit$estimate[1],rate=fit$estimate[2])
	ytest<-rgamma(ntest,shape=fit$estimate[1],rate=fit$estimate[2])
}
if (length(grep(type,"poisson"))>0){
	yfit<-dpois(as.integer(xfit),lambda=fit$estimate[1])
	ytest<-rpois(ntest,lambda=fit$estimate[1])
}
if (length(grep(type,"weibull"))>0){
	yfit<-dweibull(xfit,shape=fit$estimate[1],scale=fit$estimate[2])
}
if (length(grep(type,"cauchy"))>0)
	yfit<-dcauchy(xfit,location=fit$estimate[1],scale=fit$estimate[2])
if (length(grep(type,"chi-squared"))>0)
	yfit<-dchisq(xfit,df=fit$estimate[1],ncp=fit$estimate[2])
if (length(grep(type,"t"))>0)
	yfit<-dt(xfit,df=fit$estimate[3])
if (length(grep(type,"exponential"))>0){
	yfit<-dexp(xfit,rate=fit$estimate[1])
	ytest<-rexp(ntest,rate=fit$estimate[1])
}
if (length(grep(type,"logistic"))>0)
	yfit<-dlogis(xfit,location=fit$estimate[1],scale=fit$estimate[2])
if (length(grep(type,"geometric"))>0)
	yfit<-dgeom(xfit,prob=fit$estimate[1])
if (length(grep(type,"negative binomial"))>0){
	yfit<-dnbinom(as.integer(xfit),size=fit$estimate[1],mu=fit$estimate[2])
	ytest<-rnbinom(ntest,size=fit$estimate[1],mu=fit$estimate[2])
}
#if (length(grep(type,"lognormal"))>0)
#	yfit<-dlnorm(xfit,meanlog=fit$estimate[1],sdlog=fit$estimate[2])

kt=ks.test(x,ytest)
print(type)
print(kt)

png(outfile)
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), main=Args[6])
lines(xfit,yfit, col="red")
dev.off()

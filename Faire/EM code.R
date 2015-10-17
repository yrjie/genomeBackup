install.packages("dummies")
install.packages("ade4")
install.packages("dfoptim")

require(pscl)
require(dummies)
require(lme4)
library(ade4)
require(dfoptim)
###read the start values for the beta and b, set k as known which is estimated from fixed effect model
setwd("C:/Users/A0096342/Dropbox/r")


neg_tru_pa<-read.table(file="neg_tru_pa.csv",sep=",",header=TRUE)
neg_tru_fac<-read.table(file="neg_tru_fac.csv",sep=",",header=TRUE)
dataset<-read.table(file="Nbdata_sas.csv",sep=",",header=TRUE)
dataset$LOT<-factor(dataset$LOT)
levels(dataset$LOT)=c(1:63)
beta.in=neg_tru_pa[,1]
b.in=neg_tru_fac[,1]
Gsk.in=sd(b.in)
k=2.7313###fix k
GH=GHrule(9);absc=GH[,1];weight=GH[,2];
snp.in=c(0.7198144,0.3688318,0.673733,-0.1422731)##(w1,w2,r,mu)



###############function for the truncated negative
###density function is composed with the data set, 8 parameters beta,and the random variable for 63 lots.
tnb=function(pdata,pbeta,pb)
{
  data=pdata;beta=pbeta;b=pb#b is random variable
  beta1=beta[1];beta2=beta[2];beta3=beta[3];beta4=beta[4];
  beta5=beta[5];beta6=beta[6];beta7=beta[7];beta8=beta[8];
  Z=dummies::dummy(data$LOT);
  ita=with(data,{ita=as.vector(beta1*Width*L1+beta2*Width*L2+beta3*Width*L3+beta4*Width*L4+beta5*Width*Width*L1+
                              beta6*Width*Width*L2+beta7*Width*Width*L3+beta8*Width*Width*L4+Z%*%b)})
  lamda=exp(ita)
  t=k/(k+lamda)
  y=data$Mod
  
    value=gamma(y+k)/gamma(k)/gamma(y+1)*(t^k)*(1-t)^y/(1-t^k)
  #value gives the density variable for any dataset data, parameter beta1, and random variable,b1 
}
###density function is composed with the data set, 8 parameters beta,and the random variable for 63 lots.
logtnb=function(pdata,pbeta,pb)
{
  data=pdata;beta=pbeta;b=pb#b is random variable
  beta1=beta[1];beta2=beta[2];beta3=beta[3];beta4=beta[4];
  beta5=beta[5];beta6=beta[6];beta7=beta[7];beta8=beta[8];
  Z=dummies::dummy(data$LOT);
  ita=with(data,{ita=as.vector(beta1*Width*L1+beta2*Width*L2+beta3*Width*L3+beta4*Width*L4+beta5*Width*Width*L1+
                                 beta6*Width*Width*L2+beta7*Width*Width*L3+beta8*Width*Width*L4+Z%*%b)})
  lamda=exp(ita)
  t=k/(k+lamda)
  y=data$Mod
  
  logvalue=log(gamma(y+k))-log(gamma(k))-log(gamma(y+1))+
    k*(log(k)-log(k+lamda))+y*(log(lamda)-log(lamda+k))-log(1-t^k)
  return(logvalue)
  #value gives the density variable for any dataset data, parameter beta1, and random variable,b1 
}
######log density b
logdensb<-function(pb,psnp){
  w1=psnp[1];w2=psnp[2];r=psnp[3];mu=psnp[4]
  if((w1<=pi/2&&w1>=(-pi/2))&&(w2<=pi/2&&w2>=(-pi/2))&&r>0)
  {  
  z=(pb-mu)/r
  c1 = sin(w1); c2 = cos(w1)*sin(w2);c3=cos(w1)*cos(w2)
  a0=1.1944776*c1-0.2705981*c3; 
  a1=c2; 
  a2=-0.2705981*c1+0.6532815*c3
  pk=a0+a1*z+a2*z^2
  log.value=log((pk)^2)-(z^2)/2-log(r)
  return(log.value)}
  else
    NA
}
####density f(yi|beta) is the product likelihood for one lot

##f(yi,z)
logfyib=function(data,pbeta,pGsk,psnp,absbset)
{ w1=psnp[1];w2=psnp[2];r=psnp[3];mu=psnp[4]; 
  logdens=c(1:9)  
  for(k in 1:9){
    absb_k=absbset[k];weight_k=weight[k];absc_k=absc[k]
    logfyit_b_k=sum(logtnb(data,pbeta,absb_k))
    logden_k=logfyit_b_k+logdensb(absb_k,psnp)+log(pGsk)+log(weight_k/dnorm(absc_k))
    logdens[k]=logden_k}
  return(logdens)
  }
logbayes.b=function(data,pbeta,pb,psnp)
{ w1=psnp[1];w2=psnp[2];r=psnp[3];mu=psnp[4]; 
    logfyit_b=sum(logtnb(data,pbeta,pb))
    logdenyb=logfyit_b+logdensb(pb,psnp)
  return(-logdenyb)
}
#wit=f(yi|z)f(z)/f(yi) after change variable b to u, where u is stardard normal , 
##now f(u)is replacedby f(z)*weight
wit=function(data,pbeta,pGsk,psnp,absbset)
{ weit=c(1:9)
  logfw=logfyib(data,pbeta,pGsk,psnp,absbset)
  for(k in 1:9)
  { weit[k]=1/(sum(exp(logfw-logfw[k])))}
  return(weit)
  }
###########function for the first part of M, which is the weighted log(yi|b)
Q1=function(pbeta,absbset,rwitset)
{
  valueQ1=0
  ##begin the i step computation
  for(i in 1:63)
  {   
    data.i=get(paste0("data",i));
  logfyi.bt1=c(1:9)
       for(k in 1:9){
    absb_k=absbset[i,k]
    logfyi.bt1_k=sum(logtnb(data.i,pbeta,absb_k))
    logfyi.bt1[k]=logfyi.bt1_k
                     }
   Q1i=logfyi.bt1%*%rwitset[i,]
   valueQ1=valueQ1+Q1i;
  }
  return(-valueQ1)
}


####computations for Q2


Q2=function(psnp,absbset,rwitset)
  
{ 
  Q2=0
   for(i in 1:63)
      {   
   data.i=get(paste0("data",i)); 
  logb.t1=c(1:9)
  for(k in 1:9){
    absb_k=absb[i,k]
    weight_k=weight[k]
  logb.t1_k=logdensb(absb_k,psnp)
  logb.t1[k]=logb.t1_k}
  Q2i=logb.t1%*%rwit[i,]
   }
  Q2=Q2+Q2i
  return(-Q2)
}
  
######################the main EM function
em.NB <- function(data,snp.inits,beta.inits,b.inits,Gsk.inits,maxit=1000,tol=1e-5){
####### Initial parameter estimates
    Npdata=data
    for (i in 1:63){
    assign(paste0("data",i),subset(Npdata,Npdata$LOT==i))}
    flag <- 0
    tsnp=snp.inits;tb=b.inits;tbeta=beta.inits;tGsk=Gsk.inits;tpar=c(tbeta,tsnp)
####### Iterate between expectation and maximization steps
for(iter in 1:maxit){
   ###compute the absb needed
    absb=matrix(0,63,9)
    for(i in 1:63){
      absb[i,]=tb[i]+tGsk*absc
    };
  ###compute the weight to use
    rwit=matrix(0,63,9)
    for(i in 1:63)
    { data.i=get(paste0("data",i));
    rwit[i,]=wit(data.i,tbeta,tGsk,tsnp,absb[i,])}
  #####M step for Q1
        Q1new=optim(par=tbeta,fn=Q1,method="L-BFGS-B",
        lower=c(-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2),
        upper=c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
        absbset=absb,rwitset=rwit)
       ##get the Q1 paramters new beta
        t1beta=Q1new$par
    
  ####M step for Q2
        sw1 <- seq(-1.5, 1.5, 1.5)
        sw2 <- seq(-1.5, 1.5, 1.5)
        mb=mean(tb)
        varb=var(tb)
        gw<- expand.grid(sw1, sw2)
        lkvalue=c(1:9)
        fitr=matrix(0,9,4)
        for (g in 1:9){
        sw1=gw[[1]][g]
        sw2=gw[[2]][g]
        c1 = sin(sw1); c2 = cos(sw1)*sin(sw2);c3=cos(sw1)*cos(sw2)
        a0=1.1944776*c1-0.2705981*c3; 
        a1=c2; 
        a2=(-0.2705981*c1+0.6532815*c3)
        expect_z = 2*a0*a1+6*a1*a2;
        expect_zsq = a0^2+3*a1^2+15*a2^2+6*a0*a2;
        var_z=expect_zsq-expect_z^2
       sr=sqrt(varb/var_z)
       smu=mb-sr*expect_z   
       fit=optim(c(sw1,sw2,sr,smu),Q2,method="L-BFGS-B",lower=c(-pi/2,-pi/2,0.000001,-1),upper=c(pi/2,pi/2,1,1),absbset=absb,rwitset=rwit)
       fitr[g,]=fit$par
       lkvalue[g]=fit$value}
       min=which.min(lkvalue);
       t1snp=fitr[min,]
       ###get the Q1 new snp parameters
       t1par=c(t1beta,t1snp)
 ########## Stop iteration if the difference between the current and new estimates is less than a tolerance level
if( all(abs(tpar - t1par) < tol) ){ flag <- 1; break}
    
######## Otherwise continue iteration
else{tbeta = t1beta; tsnp=t1snp; 
#########update b
####estimate bayes estimate of bi
       b.t1=c(1:63);
       for (i in 1:63)
       { datai=get(paste0("data",i))
        fit=optimize(logbayes.b,c(-3,3),data=datai,pbeta=tbeta,psnp=tsnp)
        b.t1[i]=fit$minimum
        }
        tb=b.t1;tGsk=sd(tb)
    
     }
if(!flag) warning("Didn't converge\n")
  
    list(tbeta,tsnp,tb)
 }
}
em.NB(dataset,snp.in,beta.in,b.in,Gsk.in,maxit=1000,tol=1e-5)

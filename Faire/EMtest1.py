import numpy as np
import scipy
import scipy.optimize as op
import scipy.stats as ss

def f_sum(x, data, Z, lam):
    mu1=x[0]
    sig1=x[1]
    mu0=x[2]
    sig0=x[3]
    a=-sum(Z*scipy.log(ss.norm.pdf(data, loc=mu1, scale=sig1)*lam))
    b=-sum((1-Z)*scipy.log(ss.norm.pdf(data, loc=mu0, scale=sig0)*(1-lam)))
    ret=a+b
    #ret=-sum(Z*scipy.log(ss.norm.pdf(data, loc=mu, scale=sig)*lam))
    #print sum1,sum0
    return ret

def inv_logit(p):
    t1=scipy.exp(p)
    return t1/(1+t1)

def E_step(x1,x0,lam,data):
    a1=scipy.log(ss.norm.pdf(data, loc=x1[0], scale=x1[1])*lam)
    a0=scipy.log(ss.norm.pdf(data, loc=x0[0], scale=x0[1])*(1-lam))
    lratio=a1-a0
    return inv_logit(lratio)

eps=1e-4
data1=np.array(ss.norm.rvs(loc=20,scale=5,size=1000))
data0=np.array(ss.norm.rvs(loc=0,scale=5,size=300))
data=np.concatenate([data1,data0])
x1_old=np.array([6,1])
x0_old=np.array([-3,2])
xx_old=np.concatenate((x1_old,x0_old))
Z=np.ones(len(data))*0.5
Z[data>scipy.percentile(data,90)]=1
Z[data<scipy.percentile(data,10)]=0
#Z[1:10]=1
#Z[-10:-1]=0
lam=0.5
ans=np.zeros(len(data))
ans[0:1001]=1
diff=1
cnt=0
while diff>eps:
    lam=sum(Z)/len(data)
    res = op.minimize(f_sum,xx_old ,args=(data, Z, lam), method='Nelder-Mead', options={'xtol': 1e-8, 'disp': False})
    xx=res.x
    diff=max(abs(xx-xx_old))
    x1=xx[0:2]
    x0=xx[2:]
    Z=E_step(x1,x0,lam,data)
    xx_old=xx
    cnt+=1
print x1, x0, lam
Z[Z>0.5]=1
Z[Z<0.5]=0
print sum((ans-Z)**2),cnt

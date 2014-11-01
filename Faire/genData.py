import numpy as np
import scipy
import scipy.optimize as op
import scipy.stats as ss

def f1(x,a):
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def f2(x,a,b=10):
    return (x+100)**2

def f3(x,data):
    prob=ss.norm.pdf(data, loc=x[0], scale=x[1])
#    prob=ss.gamma.pdf(x[0],data, loc=x[1], scale=x[2])
#    prob=[scipy.exp(-((t-x[0])/x[1])**2/2)/x[1] for t in data]
    ret=-sum(scipy.log(prob))
#    if max(prob)>1:
#    	#print x, ss.gamma.pdf(x[0], 3, loc=x[1], scale=x[2])
#    	print x, prob
    return ret

def f_sum(x, data, Z, lam):
    mu=x[0]
    sig=x[1]
    ret=-sum(Z*scipy.log(ss.norm.pdf(data, loc=mu, scale=sig)*lam))
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
data1=np.array(ss.norm.rvs(loc=10,scale=5,size=1000))
data0=np.array(ss.norm.rvs(loc=0,scale=5,size=300))
data=np.concatenate([data1,data0])
print '\t'.join(['dat1']+[str(x) for x in data1])
print '\t'.join(['dat0']+[str(x) for x in data0])
#data=np.array(ss.gamma.rvs(3, loc=10,scale=5,size=3))
#print sum([ss.norm.pdf(t, loc=8.959, scale=3.28218976e-05) for t in data])
#x0=np.array([100,8])
#x0=np.array([100,8,3])
x1_old=np.array([6,1])
x0_old=np.array([-3,2])
Z=np.ones(len(data))*0.5
#Z[data>scipy.percentile(data,90)]=1
#Z[data<scipy.percentile(data,10)]=0
Z[1:10]=1
Z[-10:-1]=0
lam=0.5
ans=np.zeros(len(data))
ans[0:1001]=1
diff=1
cnt=0
while 0 and diff>eps:
    #res = op.minimize(f_sum,x1 ,args=(data, Z, lam), method='BFGS', options={'gtol': 1e-6, 'disp': True})
    #res = op.minimize(f1,[10000,521,6432,321] ,args=(data,), method='BFGS', options={'gtol': 1e-6, 'disp': True})
    res = op.minimize(f_sum,x1_old ,args=(data, Z, lam), method='Nelder-Mead', options={'xtol': 1e-8, 'disp': False})
    x1=res.x
    diff=max(abs(x1-x1_old))
    #res = op.minimize(f_sum,x0 ,args=(data, 1-Z, lam), method='BFGS', options={'gtol': 1e-8, 'disp': True})
    res = op.minimize(f_sum,x0_old ,args=(data, 1-Z, lam), method='Nelder-Mead', options={'xtol': 1e-8, 'disp': False})
    x0=res.x
    diff=max(max(abs(x0-x0_old)),diff)
    Z=E_step(x1,x0,lam,data)
    lam=sum(Z)/len(data)
    x0_old=x0
    x1_old=x1
    cnt+=1
#print x1, x0, lam
#print sum((ans-Z)**2),cnt

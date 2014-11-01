import sys,os

if len(sys.argv)<3:
    print 'Usage: eggN floorN'
    exit(1)

m=int(sys.argv[1])
n=int(sys.argv[2])
dp=[[-1]*(n+1) for i in range(m+2)]

def solve(x, y):
    if dp[x][y]>=0:
    	return dp[x][y]
    if x==1:
    	dp[x][y]=y
    if y<=1:
    	dp[x][y]=1
    if dp[x][y]>=0:
    	return dp[x][y]
    ret=-1
    for k in range(1,y):
    	tmp=max(solve(x-1, k-1), solve(x, y-k))
	if ret<0 or tmp<ret:
	    ret=tmp
    dp[x][y]=ret+1
    return dp[x][y]

def solve2(x, y):
    for i in range(1, x+1):
    	dp[i][1]=dp[i][0]=1
    for j in range(1, y+1):
    	dp[1][j]=j
    for i in range(2, x+1):
	for j in range(2, y+1):
	    now=-1
	    for k in range(1, j):
	    	tmp=max(dp[i][j-k], dp[i-1][k-1])
		if now<0 or tmp<now:
		    now=tmp
	    dp[i][j]=now+1
    return dp[x][y]

print solve2(m, n)

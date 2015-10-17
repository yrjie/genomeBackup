args=commandArgs()

set.seed(8) 

N = 1000

probs = c(.7,.85,.95)
dists = runif(N)

data = vector(length=N)
for(i in 1:N){
  if(dists[i]<probs[1]){
    data[i] = rexp(1, rate=1/60.)
  } else if(dists[i]<probs[2]){
    data[i] = rnorm(1, mean=200, sd=10)
  } else if (dists[i]<probs[3]) {
    data[i] = rnorm(1, mean=400, sd=10)
  } else {
    data[i]=rnorm(1, mean=600, sd=10)
  }
}

print(summary(data))

png("temp.png")
plot(density(data))
dev.off()


source("../Algorithm.R")
load('../../data/He/subtree.RData')

set.seed(1)
A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
m = 500
lambda.list = c(1,10,100)
p.list = c(0.1,0.2,0.3)
res =  matrix(0,nrow = length(p.list), ncol = length(lambda.list))


for(step in 1:500){
  for(i in 1:length(p.list)){
    p = p.list[i]
    for(j in 1:length(lambda.list)){
      lambda = lambda.list[j]
      X = A1[,sample(ncol(A1),m)]
      Y = runif(m,0.1,1)
      X = X[rowSums(X)!=0,]
      d = nrow(X)
      Clade = sample(d,p*d)
      X[Clade,1:(m/2)] = X[Clade,1:(m/2)]+matrix(rpois((length(Clade)*m/2), lambda),nrow = length(Clade))
      X = X*10
      X = sapply(1:m, function(i)rbinom(d, X[,i], Y[i]))
      group = c(rep(1,(m/2)), rep(2,(m/2)))
      sf.est = med(X)$sf
      res[i,j] = res[i,j]+sf_bias(sf.est,Y,group)
    }
  }
}

setwd('../table/sampling_fraction')
write.csv(res/step,'med.csv')

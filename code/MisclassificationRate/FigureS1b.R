source("../Algorithm.R")
load('../../data/He/subtree.RData')



set.seed(1)


A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>1000]
rm(A)
m = 500
eta.list = seq(0,0.1,0.01)
lambda = 2
p = 0.1
pi.1 = matrix(0,nrow = 1, ncol = length(eta.list))
error.2 = matrix(0,nrow = 1, ncol = length(eta.list))

for(step in 1:500){
  Y = c(rep(1,(m/2)),rep(2,(m/2)))
  X = A1[,sample(ncol(A1),m)]
  X = X[rowSums(X)>0,]
  Clade = order(rowSums(X),decreasing = T)[1:(p*nrow(X))]
  X[Clade,1:(m/2)] = (X[Clade,1:(m/2)])/3
  X = X*10
  X = floor(X)
  sf = runif(ncol(X))
  X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
  d = nrow(X)
  v = CStat(X)
  I0.1 = which(v>0.8)
  X0 = X[I0.1,]
  v0 = replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
  w = v[v>0.8]
  f1 = sapply(w,function(x)mean(v>x))
  f0 = sapply(w,function(x)mean(v0>x))
  pi = sum(f1*f0)/sum(f0^2)
  vord = order(v,decreasing = T)
  res = sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
  for(i in 1:length(eta.list)){
    eta = eta.list[i]
    lowerx = max(which(res[vord]<eta))
    I0 = vord[1:lowerx]
    error.2[1,i] = error.2[1,i]+mean(I0 %in% Clade)
  }
    print(step)
}
setwd('../table/accuracy')
write.csv(error.2/step,'eta_2.csv')


source("../Algorithm.R")
load('../../data/He/subtree.RData')


set.seed(1)


A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>1000]
d = nrow(A1)
rm(A)
m = 500
eta = 0.1
gamma.list = seq(0.5,0.95,0.05)
lambda = 2
p = 0.1
error.2 = matrix(0,nrow = 1, ncol = length(gamma.list))

for(step in 1:500){
  Y = runif(m,1,100)
  X = A1[,sample(ncol(A1),m)]
  X = X[rowSums(X)>0,]
  Clade = order(rowSums(X),decreasing = T)[1:(p*nrow(X))]
  X[Clade,] = t(apply(X[Clade,],1,function(x)x*Y))
  X = X*10
  X = floor(X)
  sf = runif(ncol(X))
  X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
  d = nrow(X)
  v = CStat(X)
  for(i in 1:length(gamma.list)){
    gamma = gamma.list[i]
    I0.1 = which(v>gamma)
    X0 = X[I0.1,]
    v0 = replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w = v[v>gamma]
    f1 = sapply(w,function(x)mean(v>x))
    f0 = sapply(w,function(x)mean(v0>x))
    pi = sum(f1*f0)/sum(f0^2)
    vord = order(v,decreasing = T)
    res = sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx = max(which(res[vord]<eta))
    I0 = vord[1:lowerx]
    error.2[1,i] = error.2[1,i]+mean(I0 %in% Clade)
  }
    print(step)
}
setwd('../table/accuracy')
write.csv(error.2/step,'eta_6.csv')


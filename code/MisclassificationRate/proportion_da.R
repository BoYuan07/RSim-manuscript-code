source("../Algorithm.R")
load('../../data/He/subtree.RData')

# Differential abundant taxa proportion

set.seed(1)


A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>1000]
d = nrow(A1)
#rm(A)
m = 500
prop.list = c(0.1,0.2,0.3,0.4)
lambda = 2

error = matrix(0,nrow = 1, ncol = length(prop.list))

for(step in 1:500){
    Y = c(rep(1,(m/2)),rep(2,(m/2)))
    for(i in 1:length(prop.list)){
      p = prop.list[i]
      X = A1[,sample(ncol(A1),m)]
      X = X[rowSums(X)>0,]
      Clade = sample(nrow(X),p*nrow(X))
      X[Clade,1:(m/2)] = (1+X[Clade,1:(m/2)])*lambda
      X = X*10
      sf = runif(ncol(X))
      X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
      I0.1 = rsim(X)$I0
      error[1,i] = error[1,i]+mean(I0.1 %in% Clade)
    }
    print(step)
}
setwd('../table/accuracy')
write.csv(error/step,'dag.csv')


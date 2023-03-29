source("../Algorithm.R")
load('../../data/He/subtree.RData')

# Signal Strength

set.seed(1)

A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>1000]
#rm(A)
d = nrow(A1)
m = 500
p = 0.1
lambda.list = seq(1,10,1)

error = matrix(0,nrow = 1, ncol = length(lambda.list))

for(step in 1:500){
    Clade = sample(d,p*d)
    for(i in 1:length(lambda.list)){
        lambda = lambda.list[i]
        X = A1[,sample(ncol(A1),m)]
        Y = c(rep(1,(m/2)),rep(2,(m/2)))
        X[Clade,1:(m/2)] = (X[Clade,1:(m/2)]+1)*lambda
        X = X*10
        sf = runif(ncol(X))
        X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
        I0.1 = rsim(X)$I0
        error[1,i] = error[1,i]+mean(I0.1 %in% Clade)
    }
}
 

setwd('../table/accuracy')
write.csv(error/step,'lambda.csv')

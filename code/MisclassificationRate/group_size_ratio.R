source("../Algorithm.R")
load('../../data/He/subtree.RData')

# group 1 size: group 2 size

set.seed(1)

A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>1000]
d = nrow(A1)
odds = c(0.1,0.2,0.3,0.4)
m = 500
p = 0.1
lambda = 2
error = matrix(0,nrow = 1, ncol = length(odds))

for(step in 1:500){
    Clade = sample(d,p*d)
    for(i in 1:length(odds)){
        odd = odds[i]
        X = A1[,sample(ncol(A1),m)]
        Y = c(rep(1,(m*odd)),rep(2,(m-m*odd)))
        X[Clade,1:(m*odd)] = (X[Clade,1:(m*odd)]+1)*lambda
        X = X*10
        sf = runif(ncol(X))
        X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
        I0.1 = rsim(X)$I0
        error[1,i] = error[1,i]+mean(I0.1 %in% Clade)
    }
}
setwd('../table/accuracy')
write.csv(error/step,'uneven.csv')

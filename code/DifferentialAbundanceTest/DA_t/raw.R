source("../../Algorithm.R")
load('../../../data/He/subtree.RData')

set.seed(1)
subtree = Descendants(tree, length(tree$tip.label)+351, type = "tips")[[1]]
A = A[,tree$tip.label]
drop_tips <- tree$tip.label[-subtree]
tree1 = ape::drop.tip(tree,drop_tips) %>% ape::as.phylo()
A1 = t(A[,subtree])
A1 = A1[rowSums(A1)!=0,]

d = nrow(A1)
lambda.list = c(100,200,300)
m.list = c(100,200,300)
p = 0.1

res1 =  matrix(0,nrow = length(m.list), ncol = length(lambda.list))
res2 = matrix(0,nrow = length(m.list), ncol = length(lambda.list))

for(step in 1:500){
  for(i in 1:length(m.list)){
    m = m.list[i]
    for(j in 1:length(lambda.list)){
      lambda = lambda.list[j]
      X = A1[,sample(ncol(A1),m)]
      X = X[rowSums(X)!=0,]
      Clade = sample(nrow(X),p*nrow(X))
      X[Clade,1:(m/2)] = X[Clade,1:(m/2)]+matrix(rpois((length(Clade)*m/2), lambda),nrow = length(Clade))
      X = X*lambda
      sf = exp(c(rnorm(m/2,log(1/lambda),0.01),rnorm((m-m/2),log(10/lambda),0.01)))
      sf[sf>1] = sf
      X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
      Y = c(rep(1,m/2),rep(2,m/2))
      out = t_test(X,Y,alpha = 0.1,Clade=Clade)
      res1[i,j] = res1[i,j]+out$Sensitivity
      res2[i,j] = res2[i,j]+out$FDR
    }
  }
}

setwd('../../table/t')
write.csv(res1/step,'raw_sensitivity.csv')
write.csv(res2/step,'raw_FDR.csv')

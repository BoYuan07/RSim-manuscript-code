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
lambda.list = c(50,100,150)
m.list = c(200,350,500)
p = 0.1

res1 =  matrix(0,nrow = length(m.list), ncol = length(lambda.list))
res2 = matrix(0,nrow = length(m.list), ncol = length(lambda.list))

for(step in 1:500){
  for(i in 1:length(m.list)){
    m = m.list[i]
    for(j in 1:length(lambda.list)){
      lambda = lambda.list[j]
      X = A1[,sample(ncol(A1),m)]
      X = X[rowSums(X)>0,]
      X = X[,colSums(X)>0]
      Clade = sample(nrow(X),p*nrow(X))
      Y = runif(ncol(X),1,lambda)
      X[Clade,] = sapply(1:ncol(X),function(k)X[Clade,k]+rpois(length(Clade),Y[k]))
      rand = (Y<mean(Y))
      sf = runif(m,0.5,1)*rand + runif(m,0.1,0.6)*(1-rand)
      X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
      X1 = uq(X)$P
      out = cor_test(as.matrix(X1),Y,type = "pearson",alpha = 0.1,Clade=Clade)
      res1[i,j] = res1[i,j]+out$Sensitivity
      res2[i,j] = res2[i,j]+out$FDR
    }
  }
}

setwd('../../table/pearson')
write.csv(res1/step,'uq_sensitivity.csv')
write.csv(res2/step,'uq_FDR.csv')

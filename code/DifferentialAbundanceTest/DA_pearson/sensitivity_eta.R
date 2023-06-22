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
lambda = 50
m = 200
p = 0.1
eta.list = seq(0,0.5,0.05)

res1 =  matrix(0,nrow = length(eta.list), ncol = 1)
res2 = matrix(0,nrow = length(eta.list), ncol = 1)
res3 = matrix(0,nrow = length(eta.list), ncol = 1)

for(step in 1:500){
  for(i in 1:length(eta.list)){
    eta = eta.list[i]
    X = A1[,sample(ncol(A1),m)]
    X = X[rowSums(X)>0,]
    X = X[,colSums(X)>0]
    Clade = sample(nrow(X),p*nrow(X))
    Y = runif(ncol(X),1,lambda)
    X[Clade,] = sapply(1:ncol(X),function(k)X[Clade,k]+rpois(length(Clade),Y[k]))
    rand = (Y<mean(Y))
    sf = runif(m,0.5,1)*rand + runif(m,0.1,0.6)*(1-rand)
    X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
    res = rsim(X,eta = eta)
    X1 = res$P
    out = cor_test(as.matrix(X1),Y,type = "pearson",alpha = 0.1,Clade=Clade)
    res1[i,1] = res1[i,1]+out$Sensitivity
    res2[i,1] = res2[i,1]+out$FDR
    res3[i,1] = res3[i,1]+mean(res$I0 %in% Clade)
  }
}

setwd('../../table/pearson/sensitivity')
write.csv(res1/step,'eta_sensitivity.csv')
write.csv(res2/step,'eta_FDR.csv')
write.csv(res3/step,'eta_misrate.csv')

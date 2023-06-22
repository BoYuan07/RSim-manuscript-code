## Figure 4b
## weak signal
## jobid:1-500
source("../Algorithm.R")
load('../../data/He/subtree.RData')


# weak signal

jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
simu.iter = jobid
set.seed(jobid)
A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>5000]
d = nrow(A1)
m = 500
c = 3
eta.list = seq(0,0.5,0.05)
res =  matrix(0,nrow = 1, ncol = (2*length(eta.list)))
k = 1
mode(A1) = "integer"
for(eta in eta.list){
  A1 = A1[,sample(ncol(A1),m)]
  X1 = A1[,1:floor(m/2)]
  X2 = A1[,(floor(m/2)+1):m]
  X1.1 = t(rrarefy(t(X1),colSums(X1)/c))
  ind1 = (apply(X1.1,2,function(x)sum(x!=0))<5)
  X1.1 = X1.1[,!ind1]
  X = cbind(X1.1,X2)
  X = X[rowSums(X)>0,]
  Y = c(rep(1,ncol(X1.1)),rep(2,ncol(X2)))
  X1 = rsim(X,eta=eta)$P
  res[1,k] = permanova(X1,Y)
  res[1,(k+1)] = mirkat(X1,Y)
  k = k+2
}

setwd('../table/tIe_eta')
filename = paste0("ds", "_simuiter", simu.iter, ".rds")
saveRDS(res, file = filename)

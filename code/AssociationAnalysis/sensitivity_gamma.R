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
rm(A)
m = 500
p = 0.1
gamma.list = seq(0,0.9,0.1)
res = matrix(0,nrow = 1, ncol = (length(gamma.list)*2))
eta = 0.05
k = 1
for(gamma in gamma.list){
  X = A1[,sample(ncol(A1),m)]
  X = X[rowSums(X)!=0,]
  d = nrow(X)
  Clade = order(rowSums(X),decreasing = T)[1:(p*d)]
  X[Clade,1:(m/2)] = X[Clade,1:(m/2)]*5
  sf = c(rnorm(floor(m/2),-0.7,0.05),rnorm((m-floor(m/2)),0,0.05))
  sf = exp(sf)
  sf[sf>1] = 1
  X = sapply(1:ncol(X), function(i)rbinom(nrow(X), X[,i], sf[i]))
  Y = c(rep(1,m/2),rep(2,m/2))
  X1 = rsim(X,eta=eta,gamma = gamma)$P
  res[1,k] = permanova(X1,Y)
  res[1,(k+1)] = mirkat(X1,Y)
  k = k+2
}

setwd('../table/ds_gamma')
filename = paste0("ds", "_simuiter", simu.iter, ".rds")
saveRDS(res, file = filename)

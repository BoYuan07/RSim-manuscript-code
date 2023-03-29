## Figure S5a
## strong signal
## jobid: 1-500
source("../Algorithm.R")
load('../../data/He/subtree.RData')


#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
rm(A)
m.list = c(100,300,500)
p = 0.1
res = matrix(0,nrow = 8, ncol = (length(m.list)*2))
k = 1
for(m in m.list){
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
  X.list = Normalized(X,eta=0.05,ref = (1:d)[-Clade])
  for(i in 1:length(X.list)){
      X = X.list[[i]]
      res[i,k] = permanova(X,Y)
      res[i,(k+1)] = mirkat(X,Y)
  }
  k = k+2
}

setwd('../table/ds2')
filename = paste0("ds", "_simuiter", simu.iter, ".rds")
saveRDS(res, file = filename)

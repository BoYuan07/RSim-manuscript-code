## Figure 4a
source("../Algorithm.R")
load('../../data/He/subtree.RData')

#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

A1 = t(A)
A1 = A1[rowSums(A1)!=0,]
A1 = A1[,colSums(A1)>5000]
d = nrow(A1)
m = 500
c.list = c(1,2,3)
res =  matrix(0,nrow = 8, ncol = (2*length(c.list)))
k = 1
mode(A1) = "integer"
for(c in c.list){
    A1 = A1[,sample(ncol(A1),m)]
    X1 = A1[,1:floor(m/2)]
    X2 = A1[,(floor(m/2)+1):m]
    X1.1 = t(rrarefy(t(X1),colSums(X1)/c))
    ind1 = (apply(X1.1,2,function(x)sum(x!=0))<5)
    X1.1 = X1.1[,!ind1]
    
    X = cbind(X1.1,X2)
    X = X[rowSums(X)>0,]
    Y = c(rep(1,ncol(X1.1)),rep(2,ncol(X2)))
    X.list = Normalized(X,eta=0.07)
    for(i in 1:length(X.list)){
      X = X.list[[i]]
      res[i,k] = permanova(X,Y)
      res[i,(k+1)] = mirkat(X,Y)
  }
  k = k+2
}
setwd('../table/sc')
filename = paste0("sc", "_simuiter", simu.iter, ".rds")
saveRDS(res, file = filename)

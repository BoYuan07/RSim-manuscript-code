files <- list.files(path = "./", pattern = "\\.rds$", full.names = TRUE)
data <- do.call("rbind", lapply(files, readRDS))

perm = matrix(0,nrow = 10, ncol = 3)
mirk = matrix(0,nrow = 10, ncol = 3)
for(i in 1:length(files)){
  perm = perm + (data[((i-1)*10+1):(i*10),c(1,3,5)]<0.05)
  mirk = mirk + (data[((i-1)*10+1):(i*10),c(2,4,6)]<0.05)
}
perm = perm/length(files)
mirk = mirk/length(files)
write.csv(perm, file = "../twosample/strongsig-perm.csv")
write.csv(mirk, file = "../twosample/strongsig-mirk.csv")
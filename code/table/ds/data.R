files <- list.files(path = "./", pattern = "\\.rds$", full.names = TRUE)
data <- do.call("rbind", lapply(files, readRDS))

perm = matrix(0,nrow = 8, ncol = 3)
mirk = matrix(0,nrow = 8, ncol = 3)
for(i in 1:length(files)){
  perm = perm + (data[((i-1)*8+1):(i*8),c(1,3,5)]<0.05)
  mirk = mirk + (data[((i-1)*8+1):(i*8),c(2,4,6)]<0.05)
}
perm = perm/length(files)
mirk = mirk/length(files)
write.csv(perm, file = "../twosample/ds-perm2.csv")
write.csv(mirk, file = "../twosample/ds-mirk2.csv")
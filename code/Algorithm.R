library(phangorn)
library(dplyr)
#library(ggplot2)
library(metagenomeSeq)
library(vegan)
library(MiRKAT)
library(DESeq2)
library(MicrobiomeStat)
library(matrixStats)

## Permanova ######################
permanova <- function(P, Y, p = 2) {
  if ( p == 1) {
    distP = dist(t(P), method = "minkowski", p = 1)
  } else if ( p == 2) {
    distP = dist(t(P), method = "minkowski", p = 2)
  }
  df.Y = as.data.frame(Y)
  Re = adonis2(distP~Y, data = df.Y)
  return(Re$`Pr(>F)`[1])
}

## MiRKAT #########################

mirkat <- function(P, Y, p = 2) {
  if ( p == 1) {
    distP = dist(t(P), method = "minkowski", p = 1)
  } else if ( p == 2) {
    distP = dist(t(P), method = "minkowski", p = 2)
  }
  kerlP <- D2K(as.matrix(distP))
  Re = MiRKAT(y = Y, Ks = kerlP)
  return(Re$p_values)
}


## Normalization Functions ######################################
## col: samples, row: taxa

### new method
CStat = function(X){
  d <- nrow(X)
  R <- X
  S1 <- apply(R,1,order)
  S1 <- S1 - colMeans(S1);
  S1 <- S1 / sqrt(colSums(S1^2));
  corr_s <- crossprod(S1)
  med <- as.data.frame(matrixStats::colMedians(corr_s))
  return(as.numeric(med[,1]))

}


rsim = function(X,eta=0.01,gamma = 0.8){
    d = nrow(X)
    v = CStat(X)
    I0.1 = which(v>gamma)
    X0 = X[I0.1,]
    v0 = replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w = v[v>gamma]
    f1 = sapply(w,function(x)mean(v>x))
    f0 = sapply(w,function(x)mean(v0>x))
    pi = sum(f1*f0)/sum(f0^2)
    vord = order(v,decreasing = T)
    res = sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx = max(which(res[vord]<eta))
    ref = vord[1:lowerx]
    tc.cn <- apply(X,2,function(x)sum(x[ref]))
    f.cn <- tc.cn/(mean(tc.cn))
    f.cn <- ifelse(f.cn==0,1,f.cn)
    cn.res <- scale(X,center=FALSE,scale=f.cn)
  return(list('P' = cn.res, 'I0' = ref, 'pi0'= pi, 'sf'=f.cn))
}
  
  

## true reference
true.ref <- function(X,ref){
    tc <- apply(X,2,function(x)sum(x[ref]))
    f.tc <- tc/mean(tc)
    f.tc <- ifelse(f.tc==0,1,f.tc)
    tc.res <- scale(X,center=FALSE,scale=f.tc)
    return(list('P' = tc.res,'sf'=f.tc))
}

## Upper quantile
uq<-function(X){
    dds = edgeR::calcNormFactors(as.matrix(X+1), method = "upperquartile")
    upperQ = dds*colSums(X)
    f.uq <- upperQ/mean(upperQ)
    upq.res <- scale(X,center=FALSE,scale=f.uq)
  return(list('P' = upq.res, 'sf' = f.uq))
}

## Med normalization
med<-function(X){
    logvals <- log(X)
    logvals[is.infinite(logvals)] <- NA_real_
    gm <- exp(rowMeans(logvals, na.rm=TRUE))
    med <- estimateSizeFactorsForMatrix(X, geoMeans=gm)
    f.med <- med/mean(med)
    med.res <- scale(X,center=FALSE, scale=f.med)
    return(list('P' = med.res, 'sf' = f.med))
}

## TSS
tss <- function(X){
    X <- X
    tc <- apply(X,2,sum)
    f.tc <- tc/mean(tc)
    tc.res <- scale(X,center=FALSE,scale=f.tc)
    return(list('P' = tc.res, 'sf' = f.tc))
}

## CSS
css <- function(X){
    X.css <- metagenomeSeq::newMRexperiment(X)
    X.css <- cumNorm(X.css , p=cumNormStatFast(X.css))
    P.css <- data.frame(MRcounts(X.css, norm=TRUE, log=F))
    sf <- normFactors(X.css)
    return(list('P' = P.css, 'sf' = sf))
}


## TMM
tmm <- function(X){
    f.tmm <- edgeR::calcNormFactors(as.matrix(X),method = "TMM")
    f.tmm <- f.tmm*colSums(X)
    f.tmm <- f.tmm/mean(f.tmm)
    tmm.res <- scale(X,center=FALSE,scale=f.tmm)
    return(list('P' = tmm.res, 'sf' = f.tmm))
    
}

## GMPR
GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
    # Computes the GMPR size factor
    #
    # Args:
    #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
    #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
    #   ct.min: the minimum number of counts required to calculate ratios

    #
    # Returns:
    #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
    #         NSS:   number of samples with significant sharing (> intersect.no) including itself
    
    # mask counts < ct.min
    comm[comm < ct.min] <- 0
    
    if (is.null(colnames(comm))) {
        colnames(comm) <- paste0('S', 1:ncol(comm))
    }
    
    if (trace) cat('Begin GMPR size factor calculation ...\n')
    
    comm.no <- numeric(ncol(comm))
    gmpr <- sapply(1:ncol(comm),  function(i) {
                if (i %% 50 == 0) {
                    cat(i, '\n')
                }
                x <- comm[, i]
                # Compute the pairwise ratio
                pr <- x / comm
                # Handling of the NA, NaN, Inf
                pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
                # Counting the number of non-NA, NaN, Inf
                incl.no <- colSums(!is.na(pr))
                # Calculate the median of PR
                pr.median <- colMedians(pr, na.rm=TRUE)
                # Record the number of samples used for calculating the GMPR
                comm.no[i] <<- sum(incl.no >= intersect.no)
                # Geometric mean of PR median
                if (comm.no[i] > 1) {
                    return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
                } else {
                    return(NA)
                }
            }
    )
    
    if (sum(is.na(gmpr))) {
        warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                'For these samples, their size factors are set to be NA! \n',
                'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
    }
    
    if (trace) cat('Completed!\n')
    if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
    names(gmpr) <- names(comm.no) <- colnames(comm)
    attr(gmpr, 'NSS') <- comm.no
    gmpr.res = t(t(comm) / gmpr)
    return(list('P' = gmpr.res, 'sf' = gmpr))
}

##Normalized value############################################

Normalized = function(X,eta = 0.01, seq = 5000, ref = NULL){
    rownames(X) = c(1:nrow(X))
    colnames(X) = c(1:ncol(X))
    P.cn = rsim(X,eta)$P
    P.css = css(X)$P
    P.med = med(X)$P
    P.tmm = tmm(X)$P
    P.tss = tss(X)$P
    P.uq = uq(X)$P
    P.gmpr = GMPR(X)$P
    P.ry = t(rrarefy(t(X),seq))
    res = list(P.cn,X,P.css,P.med,P.tmm,P.tss,P.uq,P.gmpr,P.ry)
    if(!is.null(ref)){
        P.true = true.ref(X,ref)$P
        res = c(res,list(P.true))
    }
    return(res)
}
est_factor = function(X,eta = 0.01, ref = NULL){
    rownames(X) = c(1:nrow(X))
    colnames(X) = c(1:ncol(X))
    sf.cn = rsim(X,eta)$sf
    sf.css = css(X)$sf
    sf.med = med(X)$sf
    sf.tmm = tmm(X)$sf
    sf.tss = tss(X)$sf
    sf.uq = uq(X)$sf
    sf.gmpr = GMPR(X)$sf
    sf.true = true.ref(X,ref)$sf
    res = data.frame('sf.est' = c(sf.cn, sf.css, sf.med, sf.tmm, sf.tss, sf.uq, sf.gmpr, sf.true),'Method' = c(rep('RSim',ncol(X)),rep('CSS',ncol(X)),rep('MED',ncol(X)),rep('TMM',ncol(X)),rep('TSS',ncol(X)),rep('UQ',ncol(X)),rep('GMPR',ncol(X)),rep('Oracle',ncol(X))))
    return(res)
}



##Differential abundance test##########################
t_test = function(X,Y,alpha,Clade,method = "BH"){
  d = nrow(X)
  data = cbind(t(X),Y)
  t_p <- sapply(1:d, function(i)t.test(data[, i] ~ Y, data = data)$p.value)
  p <- p.adjust(t_p, method=method)
  res = which(p<alpha)
  if(length(res)==0){
      fdr = 0
    }else{
      fdr = mean(!res %in% Clade)
    }
    return(list("Sensitivity" = mean(Clade %in% res), "FDR" = fdr))
}

cor_test = function(X,Y,type,alpha,Clade,method = "BH"){
  d = nrow(X)
  cor_p <- sapply(1:d, function(i)cor.test(X[i,],Y, method = type)$p.value)
  p <- p.adjust(cor_p, method=method)
  res = which(p<alpha)
  if(length(res)==0){
      fdr = 0
    }else{
      fdr = mean(!res %in% Clade)
    }
    return(list("Sensitivity" = mean(Clade %in% res), "FDR" = fdr))
}

##Sampling Fration Experiment##########################
sf_bias = function(sf.est,sf,group){
  gl = unique(group)
  index1 = which(group == gl[1])
  index2 = which(group == gl[2])
  return(mean(log((sf.est/sf)[index1]))-mean(log((sf.est/sf)[index2])))
}


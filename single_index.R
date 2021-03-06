##########################################################
############ FUNCTIONS ###################################
##########################################################
library(MCMCpack)
library(mvtnorm)

sigmod <- function(x){
  # ifelse(exp(x)==Inf,0.5,exp(x)/(1+exp(x))-0.5)
  exp(2*x)/(1+exp(2*x))-0.5
  # 0.5*x/(10+abs(x))
}

## rewrite the gpd related functions
rGpd <- function(n, loc = 0, scale = 1, shape = 0) {
  ind <- (shape == 0)
  rn <- loc + scale * (runif(n) ^ (-shape) - 1) / shape
  if (sum(ind) > 0) { rn[ind] <- loc[ind] + scale[ind] * rexp(sum(ind)) }
  return(rn)
}

dGpd <- function(x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  if (min(scale) <= 0)
    stop("invalid scale")
  d <- (x - loc) / scale
  nn <- length(d)
  scale <- rep(scale, length.out = nn)
  index <- (d >= 0 & ((1 + shape * d) > 0)) | is.na(d) 
  
  ind0 <- (shape == 0)
  rn <- log(1 / scale) - (1 / shape + 1) * log(1 + shape * d)
  if (sum(ind0,na.rm=T) > 0) { rn[ind0] <- log(1 / scale[ind0]) - d[ind0] }
  d <- ifelse(index, rn, -10000)
  d[is.na(shape)] <- -10000
  
  if (!log)
    d <- exp(d)
  d
}

pGpd <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) {
  if (min(scale) <= 0)
    stop("invalid scale")
  q <- pmax(q - loc, 0) / scale
  ind0 <- (shape == 0)
  rn <- p <- pmax(1 + shape * q, 0)
  p <- 1 - p ^ (-1 / shape)
  if (sum(ind0) > 0) { rn[ind0] <- 1 - exp(-q) }
  
  if (!lower.tail)
    p <- 1 - p
  p
}

qGpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
  if (lower.tail) 
    p <- 1 - p
  ind0 <- (shape == 0)
  rn <- loc + c(scale) * (p^(-shape) - 1)/shape
  if (sum(ind0)>0){rn[ind0] <- loc - c(scale) * log(p)} 
  
  rn
}

### truncated normal from Reich (2012)
rtnorm<-function(mn,sd=.25,a=0,b=0){
  upper<-pnorm(b,mn,sd)
  lower<-pnorm(a,mn,sd)
  if(is.matrix(mn)){
    U<-matrix(runif(prod(dim(mn)),lower,upper),dim(mn)[1],dim(mn)[2])
  }
  if(!is.matrix(mn)){
    U<-runif(length(mn),lower,upper)
  }
  return(qnorm(U,mn,sd))}

dtnorm<-function(y,mn,sd=.25,a=0,b=0){
  upper<-pnorm(b,mn,sd)
  lower<-pnorm(a,mn,sd)
  l<-dnorm(y,mn,sd,log=T)-log(upper-lower)
  return(l)}

build_b_spline = function(t, ext_knots, ind){
  w1=rep(0,length(t))
  w2=rep(0, length(t))
  b_spline <- rep(0,length(t))
  for (i in 1:length(t)){
    b_spline[i] = (ext_knots[ind-1] <= t[i])*
      (min(t[i],ext_knots[ind])-ext_knots[ind-1])
  }
  
  b_spline
}

make.eta <- function(phi){
  p <- length(phi)+1
  eta <- rep(0, p)
  sphi <- sin(phi)
  cphi <- cos(phi)
  eta[1] <- cphi[1]
  eta[p] <- prod(sphi)
  if(p>2){
    for(j in 2:(p-1)){
      eta[j] <- prod(sphi[1:(j-1)])*cphi[j]
    }
  }
  eta
}



loglike <- function(Y, Z, phi, beta1, beta2, lambda1, 
                    knots, k,
                    c1,d1,a0,s0,k.mn.prior,k.sd.prior) {
  J <- length(knots)
  n <- length(Y)
  p <- dim(Z)[2]
  
  eta <- make.eta(phi)
  X <- Z%*%eta
   
  ## basis
  lambda.fun <- matrix(0,nrow=J,ncol=n)
  for(i in 1:(J-1)){
    lambda.fun[i,] = build_b_spline(X,knots,i+1)
  }
  for(i in 1:n){
    lambda.fun[J,i] = (X[i]>knots[J])*abs(X[i]-knots[J])^k
  }
  # 
  S <- c(c(beta1) %*% rbind(1,lambda.fun))
  C <- (c(c(beta2) %*% rbind(1,lambda.fun)))
  logY <- dGpd(Y, 1, exp(S), C, log = T)
  logE <- sum(dunif(phi[1:(p-2)],0,pi,log=T)) + dunif(phi[p-1],0,2*pi,log=T)
  logB <- sum(dtnorm(beta1[-1],mn=0,sd=sqrt(rep(1/lambda1, J)),a=0,b=Inf))+
    dnorm(beta1[1],0,sd=sqrt(1/lambda1),log=T)+
    sum(dnorm(beta2,0,sd=sqrt(rep(1/lambda1, J+1)),log=T))
  logP <- dgamma(lambda1,c1,d1,log=T) + 
    dtnorm(k,k.mn.prior,sd=k.sd.prior,a=0,b=Inf)
  
  like <- sum(logY,na.rm=T) + logE + logB + logP
  
  return(like)
}


### paritial b
dlb <- function(Y, Z, phi, beta1, beta2, lambda1, knots, k) {
  
  J <- length(knots)
  n <- length(Y)
  p <- dim(Z)[2]
  
  eta <- make.eta(phi)
  X <- Z%*%eta
  
  ## basis
  lambda.fun <- matrix(0,nrow=J,ncol=n)
  for(i in 1:(J-1)){
    lambda.fun[i,] = build_b_spline(X,knots,i+1)
  }
  for(i in 1:n){
    lambda.fun[J,i] = (X[i]>knots[J])*abs(X[i]-knots[J])^k
  }
  # 
  S <- c(c(beta1) %*% rbind(1,lambda.fun))
  C <- (c(c(beta2) %*% rbind(1,lambda.fun)))
  
  dlx <- log(1+C*(Y-1)/exp(S))/C^2 - (1+1/C)*(Y-1)/(exp(S)+C*(Y-1))
  dlx[is.na(dlx)] <- 0
  
  dlb <- rbind(1,lambda.fun)%*%dlx
  dlb
}

### paritial second order b
Hb <- function(Y, Z, phi, beta1, beta2, lambda1, knots, k) {
  
  J <- length(knots)
  n <- length(Y)
  p <- dim(Z)[2]
  dlx <- rep(0,n)
  dlb <- rep(0,J)
  
  eta <- make.eta(phi)
  X <- Z%*%eta
  
  ## basis
  lambda.fun <- matrix(0,nrow=J,ncol=n)
  for(i in 1:(J-1)){
    lambda.fun[i,] = build_b_spline(X,knots,i+1)
  }
  for(i in 1:n){
    lambda.fun[J,i] = (X[i]>knots[J])*abs(X[i]-knots[J])^k
  }
  # 
  S <- c(c(beta1) %*% rbind(1,lambda.fun))
  C <- (c(c(beta2) %*% rbind(1,lambda.fun)))
  
  dl2x <- -2*log(1+C*(Y-1)/exp(S))/C^3 + 2*(Y-1)/(exp(S)+C*(Y-1))/C^2 +
    (1+1/C)*((Y-1)/(exp(S)+C*(Y-1)))^2
  dl2x[is.na(dlx)] <- 0
  
  dl2b <- rbind(1,lambda.fun)%*%diag(dl2x)%*%t(rbind(1,lambda.fun))
  dl2b
}


### paritial g
dlg <- function(Y, Z, phi, beta1, beta2, lambda1, knots, k) {
  
  J <- length(knots)
  n <- length(Y)
  p <- dim(Z)[2]
  
  eta <- make.eta(phi)
  X <- Z%*%eta
  
  ## basis
  lambda.fun <- matrix(0,nrow=J,ncol=n)
  for(i in 1:(J-1)){
    lambda.fun[i,] = build_b_spline(X,knots,i+1)
  }
  for(i in 1:n){
    lambda.fun[J,i] = (X[i]>knots[J])*abs(X[i]-knots[J])^k
  }
  # 
  S <- c(c(beta1) %*% rbind(1,lambda.fun))
  C <- (c(c(beta2) %*% rbind(1,lambda.fun)))
  
  dlx <- -1/exp(S) + (1+1/C)*((Y-1)*C/(exp(2*S)+C*(Y-1)*exp(S)))
  
  dlx[is.na(dlx)] <- 0
  
  dlb <- rbind(1,lambda.fun)%*%dlx
  dlb
}

### paritial second order
Hg <- function(Y, Z, phi, beta1, beta2, lambda1, knots, k) {
  
  J <- length(knots)
  n <- length(Y)
  p <- dim(Z)[2]
  dlx <- rep(0,n)
  dlb <- rep(0,J)
  
  eta <- make.eta(phi)
  X <- Z%*%eta
  
  ## basis
  lambda.fun <- matrix(0,nrow=J,ncol=n)
  for(i in 1:(J-1)){
    lambda.fun[i,] = build_b_spline(X,knots,i+1)
  }
  for(i in 1:n){
    lambda.fun[J,i] = (X[i]>knots[J])*abs(X[i]-knots[J])^k
  }
  # 
  S <- c(c(beta1) %*% rbind(1,lambda.fun))
  C <- (c(c(beta2) %*% rbind(1,lambda.fun)))
  
  dl2x <- 1/exp(2*S) - 
    (1+1/C)*(2*exp(S)+C*(Y-1)*C*(Y-1))/(exp(2*S)+C*(Y-1)*exp(S))^2
    
  dl2x[is.na(dlx)] <- 0
  
  dl2b <- rbind(1,lambda.fun)%*%diag(dl2x)%*%t(rbind(1,lambda.fun))
  dl2b
}

## MCMC
MCMC <- function(Y,Z,J,knots,iters=5000,burn=2500,
		 c1=J,d1=1,a0=0,s0=0.2,
		 k.mn.prior=0,k.sd.prior=0.3,
	 	 x.mn.prior=0,x.sd.prior=0.3){
  
k <- 0.1
n <- dim(Z)[1]
p <- dim(Z)[2]
J <- length(knots)
phi <- rep(1,p-1)
beta1  <- rep(0.1,J+1)
beta2  <- rep(0.01,J+1)
alpha1 <- alpha2 <- 1
lambda1 <- lambda2 <- 10
a=0;b=1;


curll <- loglike(Y, Z, phi, beta1, beta2, lambda1, 
                 knots, k,
                 c1,d1,a0,s0,k.mn.prior,k.sd.prior)

keepers <- matrix(0,iters,2) #lambda1, k
keepers.beta1 <- matrix(0,iters,J+1)
keepers.beta2 <- matrix(0,iters,J+1)
keepers.phi <- matrix(0,iters,p-1)
att <- acc <- rep(0,2)
attb1 <- accb1 <- rep(0, J+1)
attb2 <- accb2 <- rep(0, J+1)
att.phi <- acc.phi <- rep(0,p-1)
MH <- c(0.5,0.025)
MHb1 <- rep(0.01, J+1)
MHb2 <- rep(0.01, J+1)
# MHb <- 0.01
# attb2 <- accb2 <- 0
MH.phi <- rep(0.1, p-1)

for (i in 1:iters) {
  ## update lambda1
  canlambda1 = rtnorm(lambda1, MH[1],a=0,b=Inf)
  att[1] = att[1] + 1
  canll  = loglike(Y, Z, phi, beta1, beta2, canlambda1, 
                   knots, k,
                   c1,d1,a0,s0,k.mn.prior,k.sd.prior)
  R = canll - curll
  if (!is.na(exp(R))) {
    if (runif(1) < exp(R)) {
      lambda1 <- canlambda1;
      curll <- canll;
      acc[1] <- acc[1] + 1
    }
  }
  
  ## update k
  cank = rtnorm(k, MH[2],a=0,b=Inf)
  att[2] = att[2] + 1
  canll  = loglike(Y, Z, phi, beta1, beta2, lambda1, 
                   knots, cank,
                   c1,d1,a0,s0,k.mn.prior,k.sd.prior)
  R = canll - curll 
  if (!is.na(exp(R))) {
    if (runif(1) < exp(R)) {
      k <- cank;
      curll <- canll;
      acc[2] <- acc[2] + 1
    }
  }
  
  ## update eta
  canphi <- phi
  if(p>2){
    for(j in 1:(p-2)){
      att.phi[j] <- att.phi[j] + 1
      canphi[j] <- rtnorm(mn=phi[j], sd=MH.phi[j], a=0, b=pi)
      canll <- loglike(Y, Z, canphi, beta1, beta2, lambda1, 
                       knots, k,
                       c1,d1,a0,s0,k.mn.prior,k.sd.prior)
      R <- canll - curll 
      if (!is.na(exp(R))) {
        if (runif(1) < exp(R)) {
          phi[j] <- canphi[j];
          curll <- canll;
          acc.phi[j] <- acc.phi[j] + 1
        }
      }
    }
  }
  att.phi[p-1] <- att.phi[p-1] + 1
  canphi[p-1] <- rtnorm(mn=phi[p-1],sd=MH.phi[p-1],a=0,b=2*pi)
  canll <- loglike(Y, Z, canphi, beta1, beta2, lambda1, 
                   knots, k,
                   c1,d1,a0,s0,k.mn.prior,k.sd.prior)
  R <- canll - curll 
  if (!is.na(exp(R))) {
    if (runif(1) < exp(R)) {
      phi <- canphi;
      curll <- canll;
      acc.phi[p-1] <- acc.phi[p-1] + 1
    }
  }
  
  # update beta1
  canbeta1 <- beta1
  attb1[1] <- attb1[1] + 1
  canbeta1[1] <- rnorm(1,beta1[1], sd=MHb1[1])
  canll <- loglike(Y, Z, phi, canbeta1, beta2, lambda1, 
                   knots, k,
                   c1,d1,a0,s0,k.mn.prior,k.sd.prior)
  R <- canll - curll 
  if (!is.na(exp(R))) {
    if (runif(1) < exp(R)) {
      beta1[1] <- canbeta1[1];
      curll <- canll;
      accb1[1] <- accb1[1] + 1
    }
  }
  
  # canbeta1 <- beta1
  for(j in 2:(J+1)){
    attb1[j] <- attb1[j] + 1
    canbeta1[j] <- rtnorm(mn=beta1[j], sd=MHb1[j],a=0,b=Inf)
    canll <- loglike(Y, Z, phi, canbeta1, beta2, lambda1, 
                     knots, k,
                     c1,d1,a0,s0,k.mn.prior,k.sd.prior)
    R <- canll - curll 
    if (!is.na(exp(R))) {
      if (runif(1) < exp(R)) {
        beta1[j] <- canbeta1[j];
        curll <- canll;
        accb1[j] <- accb1[j] + 1
      }
    }
  }
  
  canbeta2 = beta2
  j=1
    attb2[j] <- attb2[j] + 1
    g <- dlg(Y, Z, phi, beta1, beta2, lambda1, knots, k)
    canbeta2 <- as.numeric(rmvnorm(1,mean=beta2+MHb2[j]*g,sigma=diag(J+1)*MHb2[j]))
    canll <- loglike(Y, Z, phi, beta1, canbeta2, lambda1,
                     knots, k,
                     c1,d1,a0,s0,k.mn.prior,k.sd.prior)
    R <- canll - curll
    if (!is.na(exp(R))) {
      if (runif(1) < exp(R)) {
        beta2 <- canbeta2;
        curll <- canll;
        accb2[j] <- accb2[j] + 1
      }
    }
  
  #########  TUNE THE CANDIDATE DISTRIBUTION  #######
  for (j in 1:length(acc)) {
    if (i < burn & att[j] > 50) {
      if (acc[j] / att[j] < 0.3) { MH[j] = MH[j] * 0.8 }
      if (acc[j] / att[j] > 0.6) { MH[j] = MH[j] * 1.3 }
      acc[j] = att[j] <- 0
    }
  }
  
  for (j in 1:length(acc.phi)) {
    if (i < burn & att.phi[j] > 50) {
      if (acc.phi[j] / att.phi[j] < 0.3) { MH.phi[j] = MH.phi[j] * 0.7 }
      if (acc.phi[j] / att.phi[j] > 0.6) { MH.phi[j] = MH.phi[j] * 1.3 }
      acc.phi[j] = att.phi[j] <- 0
    }
  }
  
  for (j in 1:length(accb1)) {
    if (i < burn & attb1[j] > 50) {
      if (accb1[j] / attb1[j] < 0.3) { MHb1[j] = MHb1[j] * 0.7 }
      if (accb1[j] / attb1[j] > 0.6) { MHb1[j] = MHb1[j] * 1.3 }
      accb1[j] = attb1[j] <- 0
    }
  }
  
  for (j in 1:length(accb2)) {
    if (i < burn & attb2[j] > 50) {
      if (accb2[j] / attb2[j] < 0.3) { MHb2[j] = MHb2[j] * 0.7 }
      if (accb2[j] / attb2[j] > 0.6) { MHb2[j] = MHb2[j] * 1.3 }
      accb2[j] = attb2[j] <- 0
    }
  }
  
  keepers[i,] = c(lambda1,k) #alpha1,alpha2, lambda1, lambda2, a, b, k
  keepers.phi[i,] = phi
  keepers.beta1[i,] = beta1
  keepers.beta2[i,] = beta2
    
  print(i)
}
 
  return(list(keepers=keepers, beta1=keepers.beta1, 
              phi=keepers.phi,beta2=keepers.beta2))
}
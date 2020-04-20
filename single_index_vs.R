##########################################################
############ FUNCTIONS ###################################
##########################################################
library(MCMCpack)
library(mvtnorm)

loglike <- function(Y, Z, phi, beta1, beta2, lambda1, 
                    knots, k, g,
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
  # logE <- sum(dunif(phi[1:(p-2)],0,pi,log=T)) + dunif(phi[p-1],0,2*pi,log=T)
  if(phi[p-1]<pi/4 & phi[p-1]>0){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*2/pi,10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*(dbeta(phi[p-1]*4/pi,0.01,10,log=T))/8)
  }else if(phi[p-1]<pi/2 & phi[p-1]>pi/4){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*2/pi,10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*2/pi,10,0.01,log=T)/8)
  }else if(phi[p-1]<3*pi/4 & phi[p-1]>pi/2){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*4/(3*pi),10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,0.01,10,log=T)/8)
  }else if(phi[p-1]<pi & phi[p-1]>3*pi/4){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]/pi,10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,10,0.01,log=T)/8)
  }else if(phi[p-1]<5*pi/4 & phi[p-1]>pi){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*4/(5*pi),10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,0.01,10,log=T)/8)
  }else if(phi[p-1]<3*pi/2 & phi[p-1]>5*pi/4){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*2/(3*pi),10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,10,0.01,log=T)/8)
  }else if(phi[p-1]<7*pi/4 & phi[p-1]>6*pi/4){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]*4/(7*pi),10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,0.01,10,log=T)/8)
  }else if(phi[p-1]<2*pi & phi[p-1]>7*pi/4){
    logE <- sum(g[1:(p-2)]*log(1/pi)) + 
      sum((1-g[1:(p-2)])*dbeta(phi[1:(p-2)]/(2*pi),10,0.01,log=T)) +
      sum(g[p-1]*log(1/(2*pi))) +
      sum((1-g[p-1])*dbeta(phi[p-1]*4/pi,10,0.01,log=T)/8)
  }
  
  logB <- sum(dtnorm(beta1[-1],mn=0,sd=sqrt(rep(1/lambda1, J)),a=0,b=Inf))+
    dnorm(beta1[1],0,sd=sqrt(1/lambda1),log=T)+
    sum(dnorm(beta2,0,sd=sqrt(rep(1/lambda1, J+1)),log=T))
  logP <- dgamma(lambda1,c1,d1,log=T) + 
    dtnorm(k,k.mn.prior,sd=k.sd.prior,a=0,b=Inf)+
    sum(g*log(0.3) + (1-g)*log(0.7))
  
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
  g <- rep(1,p-1)
  prq <- 0.3
  prp <- rbinom(p-1, 1, prq)
  
  
  curll <- loglike(Y, Z, phi, beta1, beta2, lambda1, 
                   knots, k,g,
                   c1,d1,a0,s0,k.mn.prior,k.sd.prior)
  
  keepers <- matrix(0,iters,11) #lambda1, k
  keepers.beta1 <- matrix(0,iters,J+1)
  keepers.beta2 <- matrix(0,iters,J+1)
  keepers.phi <- matrix(0,iters,p-1)
  att <- acc <- rep(0,2)
  attb1 <- accb1 <- rep(0, J+1)
  attb2 <- accb2 <- rep(0, J+1)
  att.phi <- acc.phi <- rep(0,p-1)
  attg <- accg <- rep(0, p-1)
  
  MH <- c(0.5,0.025)
  MHb1 <- rep(0.01, J+1)
  MHb2 <- rep(0.01, J+1)
  # MHb <- 0.01
  # attb2 <- accb2 <- 0
  MH.phi <- rep(0.1, p-1)
  MHg <- rep(0.1, p-1)
  
  
  for (i in 1:iters) {
    ## update lambda1
    canlambda1 = rtnorm(lambda1, MH[1],a=0,b=Inf)
    att[1] = att[1] + 1
    canll  = loglike(Y, Z, phi, beta1, beta2, canlambda1, 
                     knots, k, g,
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
                     knots, cank, g,
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
                         knots, k, g,
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
                     knots, k, g,
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
                     knots, k, g,
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
      canbeta1[j] <- rtnorm(mn=beta1[j], sd=MHb1[j],a=0.01,b=Inf)
      canll <- loglike(Y, Z, phi, canbeta1, beta2, lambda1, 
                       knots, k, g,
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
    gb <- dlg(Y, Z, phi, beta1, beta2, lambda1, knots, k)
    canbeta2 <- as.numeric(rmvnorm(1,mean=beta2+MHb2[j]*gb,sigma=diag(J+1)*MHb2[j]))
    canll <- loglike(Y, Z, phi, beta1, canbeta2, lambda1,
                     knots, k, g,
                     c1,d1,a0,s0,k.mn.prior,k.sd.prior)
    R <- canll - curll
    if (!is.na(exp(R))) {
      if (runif(1) < exp(R)) {
        beta2 <- canbeta2;
        curll <- canll;
        accb2[j] <- accb2[j] + 1
      }
    }
    
    # update g
    thetap <- phi
    thetap[1:(p - 2)] <- (pi / 2 - abs(phi[1:(p-2)]-pi/2))/(pi/2)
    temp <- thetap[p-1]
    temp <- temp - pi/2*(temp > pi/2) - pi/2*(temp > pi) - pi/2*(temp > 3*pi/2)
    temp <- (pi/4 - abs(temp-pi/4))/(pi/4)
    thetap[p-1] <- temp

    pr <- 0.3 / (0.3 + (1 - 0.3)*dbeta(thetap, 10, 0.01))
    
    cang = g
    for(j in 1:(p-1)){
      cang[j] = rbinom(1,1,prob=pr[j])
      canll <- loglike(Y, Z, phi, beta1, beta2, lambda1,
                       knots, k, cang,
                       c1,d1,a0,s0,k.mn.prior,k.sd.prior)
      R <- canll - curll
      if (!is.na(exp(R))) {
        if (runif(1) < exp(R)) {
          g <- cang;
          curll <- canll;
        }
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
    
    for (j in 1:length(accg)) {
      if (i < burn & attg[j] > 50) {
        if (accg[j] / attg[j] < 0.3) { MHg[j] = MHg[j] * 0.8 }
        if (accg[j] / attg[j] > 0.6) { MHg[j] = MHg[j] * 1.3 }
        accg[j] = attg[j] <- 0
      }
    }
    
    
    keepers[i,] = c(lambda1,k,g) #alpha1,alpha2, lambda1, lambda2, a, b, k
    keepers.phi[i,] = phi
    keepers.beta1[i,] = beta1
    keepers.beta2[i,] = beta2
        
    print(i)
  }
  return(list(keepers=keepers, beta1=keepers.beta1, phi=keepers.phi, beta2=keepers.beta2))
}